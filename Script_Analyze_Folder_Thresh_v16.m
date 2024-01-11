clear all;
clc

%% Main settings by the user:
path_main = '/home/min/a/salem8/lib/Fitting/Data/';
folder_name = '2023_12_13_asyn_phrodo/';
% N_measure = [11 17 28 41];  %%% Measure file number #, measure all if 0, measure a specific set if specified in a vector form [2,5,...,7]
N_measure = 11;  %%% Measure file number #, measure all if 0, measure a specific set if specified in a vector form [2,5,...,7]
Auto_threshold = 0;     %%% Threshold Slider if 0. Fixed Threshold at set value.
Show_results_only = 0;  %%% To plot the r esults only and skip the results
% filename = ["2023_11_07_neurons_trucated/well3_09";"2023_11_08_neurons_trucated/well5_04";"2023_11_08_neurons_trucated/well5_11"];

%%% Manually enter file name in to bypass search in data folder 
%%% and access the file directly, also 
%%% used with show_results_only = 1 to to to access the mat
%%% files without having the data available

%%%% Filename variable should not exist and should remain commented out if you want the code to
%%%% check the whole data folder and analyze it. Otherwise writing the folder/file name
%%%% is used when showing results only to directly go to these mat files
%%%% without having the data

%% Getting file names
warning('off','All')
path_main = strcat(path_main,folder_name);
count_file = 1;
path_ls = dir(path_main);
if (~exist('filename','var'))
    for i = 3:size(path_ls,1)
        if (path_ls(i).isdir == 1 && path_ls(i).name ~= "Properties")
            if (path_ls(i).name == "RecSettings.txt")
                continue;
            end
            filename(count_file,:) = string(path_ls(i).name);
            count_file = count_file + 1;
        end
    end
end
disp(filename)
disp(strcat('Number of files to in folder:  ',num2str(length(filename))));

if (N_measure == 0)
    N_measure  = 1:length(filename);
    sweep_count = length(filename);
else 
    sweep_count = length(N_measure);
end
disp(strcat('Number of measurements to run:  ',num2str(sweep_count)));

%% Thresholding
Thresh_array = zeros(sweep_count,1);
Mask = zeros(608,968,sweep_count);
ROI_list = cell(sweep_count,1);

if (Show_results_only == 0)
    for j = 1:sweep_count
        idx = N_measure(j);
        path = strcat(path_main,filename(idx,:));
        [segment_Obj,Threshold,Selected_ROI] = ThresholdingSlider(path,Auto_threshold);
        Thresh_array(j) = Threshold;
        Mask(:,:,j) = (segment_Obj ~= 0);
        ROI_list{j} = Selected_ROI;
    end
end

%% Lifetime fitting/viewing
[param,param_plot] = GetIniParam();
if (Show_results_only == 1)
    param.parallel = 0;
    param.SHOW_ONLY = 1;
end
if (param.parallel == 1 && param.DEBUGGING == 0)
    pool = gcp;             % Get the current parallel pool
    if ~isempty(pool)       % Parallel pool is active
        fprintf('A parallel pool is active with %d workers.\n', pool.NumWorkers);
    else
        fprintf('No parallel pool is active.\n');
        parpool('Threads')
    end
end
for i = 1:sweep_count
    if (Show_results_only == 0)
        param.THRESHOLD = Thresh_array(i);
        param.Mask = Mask(:,:,i);
        param.ROI = ROI_list{i};
    end 
    disp(strcat('Starting:  ',num2str(i),', Filename:',filename(N_measure(i),:)));
    path = strcat(path_main,filename(N_measure(i),:));
    status = SalemPixelFitting16(path, param, param_plot);
    fprintf('Done:  %d',i);
    fprintf(', Status: %d\n',status)
end

function [segment_Obj,Threshold,Selected_ROI] = ThresholdingSlider(path,Auto_threshold)
    path     = char(path);
    if (path(end) == '/')
        path = path(1:end-1);
    end
    Selected_ROI = [];
    files = dir([path, '/*.tif']);
    idx_peak = 23;  %% Manually calculated
    image_filename =  strcat(path,'/',files(idx_peak).name);
    Image_handle = Tiff(image_filename,'r');  % read given frame
    I = double(read(Image_handle)); % extract intensity data for given frame
    
    if(Auto_threshold == 0)
        [Mask,Threshold,Selected_ROI] = sliderImageDemo(I);
    else
        Threshold = Auto_threshold;
        Img = imsharpen(I,'Radius',30,'Amount',2);
        peak = max(Img(:));
        Mask = (Img > Threshold*peak);
    end

    [segment_Obj , N_Obj] = bwlabel(Mask ,8);
    %%% Remove Small objects
    segment_sizes = accumarray(segment_Obj(:)+1, 1, [N_Obj+1, 1]);
    small_segments = find(segment_sizes < 50);
    for i = 1:length(small_segments)
        mask_thresh = (segment_Obj == small_segments(i)-1);
        segment_Obj(mask_thresh) = 0;
        N_Obj = N_Obj - 1;
    end

    %%% Rename the objects
    Obj_list = unique(segment_Obj);
    for i = 2:N_Obj+1
        mask_thresh = (segment_Obj == Obj_list(i));
        segment_Obj(mask_thresh) = i-1;
    end
end
%% Plot Intensity image on top of the intensity image

function [img1,Thresh,ROI] = sliderImageDemo(Peak_Img)
    flag_ROI = 0;
    flag_sharpness = 0;
    flag_filter = 0;
    ROI = [];
    fig = uifigure;
    fig.Position = [50, 50, 1400, 500];
    fig.AutoResizeChildren = 'off';
    Thresh = 0.25;
    N_obj = 50;
    X = size(Peak_Img,1)/2;
    Y = size(Peak_Img,1)/2;
    L = 100;
    W = 100;
    A = 2;
    R = 30;
    h = -50;
    w = 150;
    w1 = 150;
    sgtitle('Thresholding/ROI','Fontsize',20,'Parent',fig);
    frame = uipanel(fig,'Position',[w+900, 20, 300, 450]);
    frame.Title = 'Controls';
    slider1 = uicontrol('Parent',fig,'Style', 'slider', 'Min', 0, 'Max', 1, 'Value', Thresh, ...
        'Position', [w1+1000, h+450, 100, 20], 'Callback', @updateImage);
    
    uicontrol('Parent',fig,'Style', 'text', 'Position', [w1+910, h+450, 90, 20], ...
        'String', 'Threshold');

    label1 = uicontrol('Parent',fig,'Style', 'text', 'Position', [w1+1100, h+450, 50, 20], ...
        'String', num2str(Thresh));

    slider_A = uicontrol('Parent',fig,'Style', 'slider', 'Min', 0, 'Max', 2, 'Value', A, ...
        'Position', [w1+1000, h+380, 100, 20], 'Callback', @updateImage);
    
    uicontrol('Parent',fig,'Style', 'text', 'Position', [w1+910, h+380, 90, 20], ...
        'String', 'A');

    label_A = uicontrol('Parent',fig,'Style', 'text', 'Position', [w1+1100, h+380, 50, 20], ...
        'String', num2str(A));

    slider_R = uicontrol('Parent',fig,'Style', 'slider', 'Min', 0, 'Max', 30, 'Value', R, ...
        'Position', [w1+1000, h+350, 100, 20], 'Callback', @updateImage);
    
    uicontrol('Parent',fig,'Style', 'text', 'Position', [w1+910, h+350, 90, 20], ...
        'String', 'R');

    label_R = uicontrol('Parent',fig,'Style', 'text', 'Position', [w1+1100, h+350, 50, 20], ...
        'String', num2str(R));

    % Create a slider
    slider_L = uicontrol('Parent',fig,'Style', 'slider', 'Min', 0, 'Max', size(Peak_Img,2), 'Value', L, ...
        'Position', [w+1000, h+200, 100, 20], 'Callback', @updateImage);
    slider_W = uicontrol('Parent',fig,'Style', 'slider', 'Min', 0, 'Max', size(Peak_Img,1), 'Value', W, ...
        'Position', [w+1000, h+230, 100, 20], 'Callback', @updateImage);
    slider_Y = uicontrol('Parent',fig,'Style', 'slider', 'Min', 0, 'Max', size(Peak_Img,1), 'Value', Y, ...
        'Position', [w+1000, h+260, 100, 20], 'Callback', @updateImage);
    slider_X = uicontrol('Parent',fig,'Style', 'slider', 'Min', 0, 'Max', size(Peak_Img,2), 'Value', X, ...
        'Position', [w+1000, h+290, 100, 20], 'Callback', @updateImage);

    uicontrol('Parent',fig,'Style', 'text', 'Position', [w+940, h+200, 50, 20], ...
        'String', 'Length');
    uicontrol('Parent',fig,'Style', 'text', 'Position', [w+940, h+230, 50, 20], ...
        'String', 'Width');
    uicontrol('Parent',fig,'Style', 'text', 'Position', [w+940, h+260, 50, 20], ...
        'String', 'Y');
    uicontrol('Parent',fig,'Style', 'text', 'Position', [w+940, h+290, 50, 20], ...
        'String', 'X');

    label_L = uicontrol('Parent',fig,'Style', 'text', 'Position', [w+1100, h+200, 50, 20], ...
        'String', num2str(L));
    label_W = uicontrol('Parent',fig,'Style', 'text', 'Position', [w+1100, h+230, 50, 20], ...
        'String', num2str(W));
    label_Y = uicontrol('Parent',fig,'Style', 'text', 'Position', [w+1100, h+260, 50, 20], ...
        'String', num2str(Y));
    label_X = uicontrol('Parent',fig,'Style', 'text', 'Position', [w+1100, h+290, 50, 20], ...
        'String', num2str(X));


    %%% Objects filter
    slider_N = uicontrol('Parent',fig,'Style', 'slider', 'Min', 0, 'Max', 150, 'Value', N_obj, ...
        'Position', [w+1000, h+130, 100, 20], 'Callback', @updateImage);
    
    uicontrol('Parent',fig,'Style', 'text', 'Position', [w1+910, h+130, 90, 20], ...
        'String', 'Obj. Fitler');

    label_N = uicontrol('Parent',fig,'Style', 'text', 'Position', [w1+1100, h+130, 50, 20], ...
        'String', num2str(N_obj));

    %%% Enable initially
    slider_L.Enable = 'off';
    slider_W.Enable = 'off';
    slider_X.Enable = 'off';
    slider_Y.Enable = 'off';
    slider_N.Enable = 'off';
    slider_A.Enable = 'on';
    slider_R.Enable = 'on';

    checkbox = uicheckbox(fig, 'Text', 'Custom ROI');
    checkbox.Position = [w+950, h+320, 100, 22];
    checkbox.ValueChangedFcn = @toggleSliderVisibility;
    uicontrol('Parent',fig,'Style', 'pushbutton','String',"Start Analyzing",'Position', [w+1070, h+80, 100, 30],'Callback',@SetT)

    checkbox2 = uicheckbox(fig, 'Text', 'Sharpness Mask');
    checkbox2.Position = [w+950, h+410, 100, 22];
    checkbox2.ValueChangedFcn = @toggleSharpness;

    checkbox3 = uicheckbox(fig, 'Text', 'Objects Filter');
    checkbox3.Position = [w+950, h+160, 100, 22];
    checkbox3.ValueChangedFcn = @toggleFilter;

    img0 = Peak_Img;
    img1 = Thresh_Image(Thresh,Peak_Img);
    bx1 = axes('Parent', fig);
    bx2 = axes('Parent', fig);
    bx1.Position = [0.05 0.25 0.3 0.5];
    bx2.Position = [0.4 0.25 0.3 0.5];
    imgHandle0 = imagesc(bx1,img0);
    imgHandle1 = imagesc(bx2,img1);
    
    rectangleHandle1 = rectangle(bx1,'Position', [X, Y, W, L], 'EdgeColor', 'b', 'FaceColor', 'none','Visible','off');
    rectangleHandle2 = rectangle(bx2,'Position', [X, Y, W, L], 'EdgeColor', 'b', 'FaceColor', 'none','Visible','off');

    colormap(bx1 , gray(256));
    colormap(bx2 , gray(2));
    colorbar(bx1);
    title('Peak Intensity','Parent',bx1);
    title('Thresholded Image','Parent',bx2);
    cb = colorbar(bx2);
    clim(bx2,[0,1]);
    cb.Ticks = [1,2] ; 
    cb.YTick = cb.Limits;
    cb.TickLabels = {'Off','On'} ;
    
    uiwait(fig)

function toggleSliderVisibility(~,~)
    if checkbox.Value
        slider_L.Enable  = 'on';
        slider_W.Enable  = 'on';
        slider_X.Enable  = 'on';
        slider_Y.Enable  = 'on';

        flag_ROI = 1;
        L = get(slider_L, 'Value');
        W = get(slider_W, 'Value');
        X = get(slider_X, 'Value');
        Y = get(slider_Y, 'Value');        
        rectangleHandle1.Visible = 'on';
        rectangleHandle2.Visible = 'on';
    else
        slider_L.Enable  = 'off';
        slider_W.Enable  = 'off';
        slider_X.Enable  = 'off';
        slider_Y.Enable  = 'off';

        flag_ROI = 0;
        rectangleHandle1.Visible = 'off';
        rectangleHandle2.Visible = 'off';
    end
end


function toggleSharpness(~,~)
        if checkbox2.Value
            slider_A.Enable  = 'on';
            slider_R.Enable  = 'on';
            flag_sharpness = 1;
            A = get(slider_A, 'Value');
            R = get(slider_R, 'Value');
        else
            flag_sharpness = 0;
            slider_A.Enable  = 'off';
            slider_R.Enable  = 'off';
            set(imgHandle0, 'CData', Peak_Img);
        end
        updateImage
    end

function toggleFilter(~,~)
    if checkbox3.Value
        slider_N.Enable  = 'on';
        flag_filter  = 1;
        N = get(slider_N, 'Value');
        set(label_N,'String',num2str(N));
    else
        flag_filter  = 0;
        slider_N.Enable  = 'off';
        set(imgHandle0, 'CData', Peak_Img);
    end
    updateImage
end

    function updateImage(~,~)
        Thresh = get(slider1, 'Value');    
        img1 = Thresh_Image(Thresh,Peak_Img);
        if (flag_sharpness == 1)
            img0 = Get_Peak_Image;
            img0 = img0/max(img0(:));
            set(imgHandle0, 'CData', img0);
        end
        if (flag_filter == 1)
            N_obj = get(slider_N, 'Value');
            set(label_N,'String',num2str(N_obj));
        end

        L = round(get(slider_L, 'Value'));
        W = round(get(slider_W, 'Value'));
        X = round(get(slider_X, 'Value'));
        Y = round(get(slider_Y, 'Value'));
        set(imgHandle1, 'CData', img1);
        if (flag_ROI == 1)
            rectangleHandle1.Position = [X,Y,W,L];
            rectangleHandle2.Position = [X,Y,W,L];
        end
        set(label1, 'String', Thresh);
        set(label_L, 'String', L);
        set(label_W, 'String', W);
        set(label_X, 'String', X);
        set(label_Y, 'String', Y);
        set(label_A, 'String', A);
        set(label_R, 'String', R);
    end

    function SetT(~,~)
        Thresh = get(slider1, 'Value');    
        img1 = Thresh_Image(Thresh,Peak_Img);
        if (flag_ROI == 1)
            L = round(get(slider_L, 'Value'));
            W = round(get(slider_W, 'Value'));
            X = round(get(slider_X, 'Value'));
            Y = round(get(slider_Y, 'Value'));
            ROI = [X,Y,W,L];
        else
            ROI = [];
        end
        delete(fig)
    end

    function Mask = Thresh_Image(T,Peak_Img)    
        Image = Peak_Img;
        if (flag_sharpness == 1)
            A = get(slider_A, 'Value');
            R = get(slider_R, 'Value');
            Image = imsharpen(Image,'Radius',R,'Amount',A);
        end
        peak = max(Image(:));
        Mask = (Image > T*peak);
        if (flag_filter == 1)
            N_filter = round(get(slider_N, 'Value'));
            Mask  = bwareaopen(Mask ,N_filter,8);
        end
    end
    
    function Img = Get_Peak_Image(~,~)    
        A = get(slider_A, 'Value');
        R = get(slider_R, 'Value');
        Img = imsharpen(Peak_Img,'Radius',R,'Amount',A);
    end
end