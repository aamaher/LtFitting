clear all
close all
%% Main settings by the user:
path_main = '/home/min/a/salem8/lib/Fitting/Data/';
folder_name = ["2023_12_12_asyn_phrodo_plate1";"2023_12_13_asyn_phrodo_plate1";"2023_12_14_asyn_phrodo_plate1";"2023_12_15_asyn_phrodo_plate1"];
% folder_name = ["2023_12_13_asyn_phrodo_plate1";"2023_12_14_asyn_phrodo_plate1";"2023_12_15_asyn_phrodo_plate1"];
% folder_name = ["2023_12_13_asyn_phrodo_plate1";"2023_12_14_asyn_phrodo_plate1"];

N_measure = 11;         %%% Measure file number #, measure all if 0
Mode = 1;   %%% 1) Get Mask, 2) Lifetime Analysis, 3) Plot results 4) Re-annotate 5) Interactive

%% Finding File names
N_days = length(folder_name);
warning('off','All')
for idx_Day = 1:N_days      %%% Find file names inside each folder
    path_current = strcat(path_main,folder_name(idx_Day));
    path_ls = dir(path_current);
    idx_Measurement = 1;
    for i = 3:size(path_ls,1)
        if (path_ls(i).isdir == 1)
            if (path_ls(i).name == "RecSettings.txt" || path_ls(i).name == "Properties")
                continue;
            end    

            filename{idx_Day,idx_Measurement} = string(path_ls(i).name);
            idx_Measurement = idx_Measurement+1;
        end
    end           
end  
if (N_measure == 0)
    N_measure  = 1:size(filename,2);
    sweep_count = length(filename);
else 
    sweep_count = length(N_measure);
end

%% Segmenting Object based on thresholding, getting the masks

if (Mode  == 1 || Mode == 4)
    Mask_file = cell(N_days,1);
    path = cell(N_days,1);
    for idx_Measurement = 1:sweep_count
        for idx_Day = 1:N_days
            path{idx_Day} = strcat(path_main,folder_name(idx_Day,1),'/',filename{idx_Day,N_measure(idx_Measurement)});
            Mask_file{idx_Day} = strcat('Results/Masks/Mask_',folder_name(idx_Day,1),'_',filename{idx_Day,N_measure(idx_Measurement)},'.mat');
        end
        if (Mode == 1)
            [Label_Mask_cells,Type_Mask_cells,N_Obj,Img_Big_Seg,Img_Big_Full] = Obj_segmentation(path,1,[]);
            for idx_Day = 1:N_days
                Label_Mask = Label_Mask_cells{idx_Day};
                Type_Mask = Type_Mask_cells{idx_Day};
                save(Mask_file{idx_Day},'Label_Mask','Type_Mask','N_Obj','Img_Big_Seg','Img_Big_Full','disp_mask','-mat');
            end
        elseif (Mode == 4)
            [Label_Mask_cells,Type_Mask_cells,N_Obj,Img_Big_Seg,Img_Big_Full] = Obj_segmentation(path,2,[]);
            for idx_Day = 1:N_days
                Label_Mask = Label_Mask_cells{idx_Day};
                Type_Mask = Type_Mask_cells{idx_Day};
                save(Mask_file{idx_Day},'Label_Mask','Type_Mask','N_Obj','Img_Big_Seg','Img_Big_Full','disp_mask','-mat');
            end
        end
        
    end
end

%% Sweep across the objects and Perform the lifetime fitting for each Object.

if (Mode == 2)
    [param,param_plot] = GetIniParam();
    param_plot.plot_flag = 0;      %% No need to plot during the file, we do this later anyway
    if (param.parallel == 1 && param.DEBUGGING == 0)
        pool = gcp;             % Get the current parallel pool
        if ~isempty(pool)       % Parallel pool is active
            fprintf('A parallel pool is active with %d workers.\n', pool.NumWorkers);
        else
            fprintf('No parallel pool is active.\n');
            parpool('Threads')
        end
    end
    disp(['Total Objects: ' num2str(sweep_count) ', Measured in: ' num2str(N_days) ' days, Total Measurements' num2str(sweep_count*sweep_count)])
    for idx_Measurement = 1:sweep_count
        fprintf('\nStarting File: %d',idx_Measurement);
        for idx_Day = 1:N_days
            fprintf('\nStarting Day: %d',idx_Measurement);
            path = strcat(path_main,folder_name(idx_Day,1),'/',filename{idx_Day,N_measure(idx_Measurement)});
            disp(strcat('Starting:  ',num2str(idx_Day),', Filename:',filename{idx_Day,N_measure(idx_Measurement)}));
            % load(strcat('Results/Masks/Mask_',folder_name(idx_Day,1),'_',filename{idx_Day,N_measure(idx_Measurement)},'.mat'),'Label_Mask','N_Obj','-mat');
            load(strcat('Results/Masks/Mask_',folder_name(idx_Day,1),'_',filename{idx_Day,N_measure(idx_Measurement)},'.mat'),'Label_Mask','-mat');
            param.Mask = Label_Mask;
            status = SalemPixelFitting16(path, param,param_plot);
            fprintf('Done day:  %d\%d, Status: %d\n',idx_Day,N_days,status);
        end
        fprintf('Day: %d finished!\n',idx_Measurement);
    end
end
%% Plot the Results
if (Mode == 3)

    if (param.cont == 1)
        str_cont = 'C';
    else
        str_cont = 'D';
    end

    for idx_Measurement = 1:sweep_count
        tau_map_all  = cell(N_days,1);
        chi_map_all  = cell(N_days,1);
        Peak_Img_all = cell(N_days,1);
        Mask         = cell(N_days,1);
        for idx_Day = 1:N_days
            path = strcat(path_main,filename{idx_Day,N_measure(idx_Measurement)});
            if (path(end) == '/')
                path = path(1:end-1);
            end
            
            filename_saving_temp = split(strip((path)),'/');
            filename_saving = strcat(string(filename_saving_temp(end-1)),string(filename_saving_temp(end)));
            filename_saving_results = strcat('Results/RESULTS_EXP',num2str(param.order),str_cont,'_',filename_saving,'_M',num2str(param.Method),'.mat');
            
            % load(strcat('Results/Masks/Mask_',folder_name(idx_Day,1),'_',filename{idx_Day,N_measure(idx_Measurement)},'.mat'),'Label_Mask','N_Obj','-mat');
            load(strcat('Results/Masks/Mask_',folder_name(idx_Day,1),'_',filename{idx_Day,N_measure(idx_Measurement)},'.mat'),'Label_Mask','-mat');
            load(filename_saving_results,'I_avg','Peak_Img','tau_map','chi_map','HorizontalPixels','VerticalPixels')
            tau_map_all {idx_Day} = tau_map;
            chi_map_all {idx_Day} = chi_map; 
            Peak_Img_all{idx_Day} = Peak_Img;
            Mask{idx_Day} = Label_Mask;
        end
        plotObjDays4;
    end
    
end

%% Interactive Mode
if (Mode == 5)
    for idx_Measurement = 1:sweep_count
        plotObjLabels5_interactive(N_days,path_main,folder_name,filename,N_measure,idx_Measurement)
    end
end


function [Final_Label_Mask,Final_Type_Mask,N_Obj,Img_Big_Seg_temp,Img_Big_Full,disp_mask] = Obj_segmentation(path,Mode,file_name)
    N_days = size(path,1);
    idx_peak = 33;
    I_Peak = zeros(1216/2,1936/2,N_days);

    for idx_day = 1:N_days
        path_current   = char(path{idx_day,1});
        if (path_current(end) == '/')
            path_current = path_current(1:end-1);
        end
        files = dir([path_current, '/*.tif']);
        image_filename =  strcat(path_current,'/',files(idx_peak).name);
        Image_handle = Tiff(image_filename,'r');  % read given frame
        Im = read(Image_handle); % extract intensity data for given frame
        I_Peak(:,:,idx_day) = Im;
    end

    Img_temp = I_Peak(:,:,2);
    Img_trunc = (Img_temp<8);
    Img_trunc = bwmorph(Img_trunc,'bridge',inf);
    Img_trunc = bwareaopen(Img_trunc,50,4);
    Img_trunc = 1-Img_trunc;
    Img_trunc = bwareaopen(Img_trunc,50,4);
    Final_Mask = cell(N_days,1);
    Final_Mask_Annotated = cell(N_days,1);
    
    %% Thresholding, Registeration, Masks mixing, and Segmentation    
    if (Mode == 1)
        
        [Mask,disp_mask] = sliderImageMask(I_Peak,Img_trunc);
        d_max = max(cell2mat(disp_mask));
        d_min = min(cell2mat(disp_mask));
        x_new = size(I_Peak,1)+d_max(2)-d_min(2);
        y_new = size(I_Peak,2)+d_max(1)-d_min(1);
        Mask_Big_Full = zeros(x_new,y_new);
        Img_Big_temp = zeros(x_new,y_new);
        Img_Big_Full = zeros(x_new,y_new);
        for i = 1:N_days
            disp_xy = -disp_mask{i} + [d_max(1) d_max(2)];
            Img_Big_temp(1+disp_xy(2):disp_xy(2)+608, 1+disp_xy(1):disp_xy(1)+968) = Img_trunc.*I_Peak(:,:,i);
            Img_Big_Full = Img_Big_Full + Img_Big_temp;
            Img_Big_temp(1+disp_xy(2):disp_xy(2)+608, 1+disp_xy(1):disp_xy(1)+968) = Mask(:,:,i);
            Mask_Big_Full = Mask_Big_Full|Img_Big_temp;
        end
        Img_Big_Full = Img_Big_Full/max(Img_Big_Full(:));
    
    %%%  Using Watershed Segmentation Method
        bw  = Mask_Big_Full;
        D = -bwdist(~bw,'quasi-euclidean');
        Ld = watershed(D);
        bw2 = bw;
        bw2(Ld == 0) = 0;
        mask_em = imextendedmin(D,2);
        D2 = imimposemin(D,mask_em);
        Ld2 = watershed(D2);
        bw3 = bw;
        bw3(Ld2 == 0) = 0;
        Img_Big_Seg_temp = ~bwareaopen(~bw3, 40);
    elseif (Mode == 2)
        load(file_name,'Img_Big_Seg_temp','Img_Big_Full','disp_mask','-mat')
    end

%%  Annotate
    [Label_Mask,Type_Mask,N_Obj] = sliderAnnotate(Img_Big_Seg_temp,Img_Big_Full);

%%  Spread the segmented mask to all images
    for i = 1:N_days
        disp_xy = -disp_mask{i} + [d_max(1) d_max(2)];
        Final_Label_Mask{i} = Img_trunc.*Label_Mask(1+disp_xy(2):disp_xy(2)+608, 1+disp_xy(1):disp_xy(1)+968);
        Final_Type_Mask{i} = Img_trunc.*Type_Mask(1+disp_xy(2):disp_xy(2)+608, 1+disp_xy(1):disp_xy(1)+968);
    end

%% Show the results of the registeration. MAKE SURE TO TRUNCATE THE
    gcf = figure;
    subplots_pairs = cell(N_days,1);
    subplot_ptr = cell(N_days,1);
    for i= 1:N_days
        Img_temp = I_Peak(:,:,i);
        Img_overlay = labeloverlay(Img_temp/max(Img_temp(:)),Img_trunc.*Final_Label_Mask{i},'Colormap',jet(255));
        subplots_pairs{i} = subplot(2, 2, i,'Parent',gcf);
        
        subplot_ptr{i} = imshow(Img_overlay,'Parent',subplots_pairs{i});
        % clim(subplots_pairs{i},[1,N_days])
    end

%% Plot Intensity image on top of the intensity image
    function [Obj_Mask,disp_mask] = sliderImageMask(I,Img_trunc)
        fig = uifigure;
        fig2 = uifigure;
        gcf = figure(7);
        ax{1}     = axes(gcf);
        fig.Position = [50, 400, 1500, 350];
        fig.AutoResizeChildren = 'off';
        N_days          = size(I,3);
        sliders         = cell(N_days,1);
        sliders_filter  = cell(N_days,1);
        label           = cell(N_days,1);
        label_filter    = cell(N_days,1);
        T               = cell(N_days,1);
        subplots_pairs  = cell(N_days,1);
        imgHandle0      = cell(N_days,1);
        Used_mask       = cell(N_days,1);
        Used_mask_temp  = cell(N_days,1);
        imgHandleR      = cell(N_days,1);
        sliders_x       = cell(N_days,1);
        sliders_y       = cell(N_days,1);
        checkbox        = cell(N_days,1);
        rb              = cell(N_days,1);
        X_pos           = zeros(N_days,1);
        Y_pos           = zeros(N_days,1);
        Thresh = 0.3;
        N_filter = 30;
        h = -20;
        w =  50;
        h2 = -10;
        w2 = 100;
        x0=10;
        y0=200;
        width=960;
        height=600;
        fig2.Position = [20, 50, 1500, 250];
        fig2.AutoResizeChildren = 'off';
        
        uicontrol('Parent',fig,'Style', 'text', 'Position', [w, h+60, 100, 20],'String', 'Threshold');
        uicontrol('Parent',fig,'Style', 'text', 'Position', [w, h+30, 100, 20],'String', 'Filter');
        set(gcf,'Color','k')
        set(gcf,'position',[x0,y0,width,height])

        bg = uibuttongroup(fig2,'Position',[20 80 123 130],'Title','Reference');
        for i = 1:N_days    
            sliders{i}        = uicontrol('Parent',fig,'Style', 'slider', 'Min', 0, 'Max', 1, 'Value', Thresh  ,'Position', [w+100+(i-1)*360, h+60, 110, 20] , 'Callback', @updateImage);           
            sliders_filter{i} = uicontrol('Parent',fig,'Style', 'slider', 'Min', 0, 'Max', 120, 'Value', N_filter,'Position', [w+100+(i-1)*360, h+30, 110, 20] , 'Callback', @updateImage);
            label{i} = uicontrol('Parent',fig,'Style', 'text', 'Position', [w+220+(i-1)*360, h+60, 60, 20],'String', num2str(Thresh));
            label_filter{i} = uicontrol('Parent',fig,'Style', 'text', 'Position', [w+220+(i-1)*360, h+30, 60, 20],'String', num2str(N_filter));

            T{i} = Thresh;      % Variable used for later
            N_f{i} = N_filter;
            img0 = I(:,:,i);
            img0 = img0/max(img0(:));
            ax{i}.Position = [0.1300 0.1100 0.7750 0.8150];
            I(:,:,i) = img0;
            Mask{i} = Thresh_Image(Thresh,I(:,:,i),N_filter);
            Used_mask{i} = Mask{i};
            Used_mask_temp{i} = Mask{i};
            Img_overlay = labeloverlay(img0,i*Mask{i},'Colormap',jet(N_days));
            ax2{i}     = axes(fig);
            imgHandle0{i} = imshow(Img_overlay,'Parent',ax2{i});
            ax2{i}.Position = [0.25*(i-1) 0.25 0.3 1];
            xlim(ax2{i},[90,770]);
            axis(ax2{i},'square')
            set(ax2{i},'XColor', 'none','YColor','none')
            colormap(ax2{i}, gray(256));
            title(['Day:' num2str(i)] ,'Parent',ax2{i});
            if  (i~=1)
                ax{i} = copyobj(ax{1},gcf);
                ax{i}.UserData = linkprop([ax{1},ax{i}],{'x','y','Position'});
            end
            set(ax{i},'XColor', 'none','YColor','none')
            img_temp = Thresh_Image(Thresh,I(:,:,i),N_filter);
            imgHandleR{i} = imagesc(ax{i},i*img_temp, 'AlphaData', img_temp);
            clim(ax{i},[1,N_days]);
            ax{i}.Visible = 'off';
            colormap(ax{i},jet(N_days));
            sliders_x{i} = uislider('Parent',fig2, 'Limits', [-400 400], 'Value', 0, 'Position', [w2+110+(i-1)*310, h2+60, 150, 20] , 'ValueChangingFcn', @updateImageReg,'ValueChangedFcn', @updateImageReg);
            sliders_y{i} = uislider('Parent',fig2, 'Orientation','Vertical','Limits', [-400 400], 'Value', 0, 'Position', [w2+100+(i-1)*310, h2+80, 20, 150] , 'ValueChangingFcn',  @updateImageReg,'ValueChangedFcn', @updateImageReg);
            checkbox{i} = uicheckbox(fig2, 'Text', 'Show Mask','Position',[w2+180+(i-1)*310, h2+150,  100, 20],'ValueChangedFcn',@toggleMask,'Value',1);
            rb{i} = uiradiobutton(bg,'Position',[10, 100-20*i,  100, 20],'Text',strcat('Day ',num2str(i)));
            colorbar(ax{i},'Ticks',1:N_days,'Color','White','FontSize',17,'TickLabels',string(split(num2str(1:N_days))));
        end
        sliders_x{1}.Enable = 0;
        sliders_y{1}.Enable = 0;

        uicontrol('Parent',fig,'Style', 'pushbutton','String',"Update Reg",'Position', [w+350+(i-1)*360 h+55 100 35],'Callback',@SetReg)
        uicontrol('Parent',fig2,'Style', 'pushbutton','String',"Start Analyzing",'Position', [w+350+(i-1)*310 70 100 35],'Callback',@SetT)
        uicontrol('Parent',fig2,'Style', 'pushbutton','String',"Auto Register"  ,'Position', [w+350+(i-1)*310 120 100 35],'Callback',@SetRegAuto)
        uicontrol('Parent',fig2,'Style', 'pushbutton','String',"Change Reference"  ,'Position', [30 15 100 35],'Callback',@SetRef)
        uiwait(fig)

    %% Functions
    function SetRegAuto(~,~)
        %%% Find Reference
        for i = 1:N_days
            if rb{i}.Value
                Ref = i;
            end
        end

        for i = 1:N_days
            if ~rb{i}.Value
                tform = imregcorr(double(Used_mask_temp{i}),double(Used_mask_temp{Ref}),"translation");
                x_axis = tform.Translation(1);
                y_axis = tform.Translation(2);
                Used_mask_temp{i} = logical(imtranslate(Used_mask{i},[x_axis y_axis]));
                set(imgHandleR{i}, 'CData'    , i*Used_mask_temp{i});
                set(imgHandleR{i}, 'AlphaData', Used_mask_temp{i});
                sliders_x{i}.Value = x_axis;
                sliders_y{i}.Value = y_axis;
            end
        end
        
    end

    function SetRef(~,~)
        for i = 1:N_days
            if rb{i}.Value
                Used_mask_temp{i} = Used_mask{i};
                sliders_x{i}.Enable = 0;
                sliders_y{i}.Enable = 0;
                sliders_x{i}.Value = 0;
                sliders_y{i}.Value = 0;
                if checkbox{i}.Value
                    set(imgHandleR{i}, 'AlphaData', Used_mask_temp{i});
                    set(imgHandleR{i}, 'CData', i*Used_mask_temp{i});
                end
            else
                sliders_x{i}.Enable = 1;
                sliders_y{i}.Enable = 1;
            end
        end
    end   
    
    function toggleMask(~,~)
        for i = 1:N_days 
            if checkbox{i}.Value
                set(imgHandleR{i}, 'AlphaData', Used_mask_temp{i});
                set(imgHandleR{i}, 'CData', i*Used_mask_temp{i});
            else
                set(imgHandleR{i}, 'AlphaData', 0);
            end
        end
    end

    function updateImage(~,~)
        for i = 1:N_days
            Thresh = get(sliders{i}, 'Value');
            N_filter = round(get(sliders_filter{i}, 'Value'));
            if (T{i} ~= Thresh || N_f{i} ~= N_filter)
                T{i} = Thresh;
                N_f{i} = N_filter;
                Mask_temp = Thresh_Image(Thresh,I(:,:,i),N_filter);
                Img_overlay = labeloverlay(I(:,:,i),i*Mask_temp,'Colormap',jet(N_days));
                set(imgHandle0{i}, 'CData', Img_overlay);
                set(label{i}, 'String', Thresh);
                set(label_filter{i}, 'String', N_filter);
            end
        end
    end
 
    function updateImageReg(~,~)
        for i = 1:N_days
            x_axis = round(get(sliders_x{i}, 'Value'));
            y_axis = round(get(sliders_y{i}, 'Value'));
            if (X_pos(i) ~= x_axis || Y_pos(i) ~= y_axis)
                X_pos(i) = x_axis;
                Y_pos(i) = y_axis;
                Used_mask_temp{i} = imtranslate(Used_mask{i},[x_axis y_axis]);
                set(imgHandleR{i}, 'CData'    , i*Used_mask_temp{i});
                set(imgHandleR{i}, 'AlphaData', Used_mask_temp{i});
                colormap(ax{i},jet(N_days));
            end
        end
    end

    function SetReg(~,~)
        for i = 1:N_days
            Thresh   = get(sliders{i}, 'Value');
            N_filter = round(get(sliders_filter{i}, 'Value'));
            Used_mask{i} = Thresh_Image(Thresh,I(:,:,i),N_filter);
            set(imgHandleR{i}, 'CData', i*Used_mask{i});
            set(imgHandleR{i}, 'AlphaData', Used_mask{i});
        end
    end

    function SetT(~,~)
        disp_mask = cell(N_days,1);
        for i = 1:N_days
            Thresh = get(sliders{i}, 'Value');
            N_filter = round(get(sliders_filter{i}, 'Value'));
            Obj_Mask(:,:,i) = Thresh_Image(Thresh,I(:,:,i),N_filter);
            if rb{i}.Value
                disp_mask{i} = [0 0];
            else
                x_axis = round(get(sliders_x{i}, 'Value'));
                y_axis = round(get(sliders_y{i}, 'Value'));
                disp_mask{i} = [-x_axis -y_axis];
            end
        end
        
        delete(fig)
        delete(fig2)
        close(gcf)
    end
    
    function mask_temp = Thresh_Image(T,Peak_Img,N_filter)
        Peak_Img  = imsharpen(Peak_Img,'Radius',30,'Amount',2);
        Peak_Img  = adapthisteq(Img_trunc.*Peak_Img/max(Peak_Img(:))); 
        mask_temp = imbinarize(Peak_Img,T);
        mask_temp = bwareaopen(mask_temp ,N_filter,8);
        mask_temp = bwmorph(mask_temp ,'bridge',inf);
        mask_temp = bwmorph(mask_temp ,'fill',1);       
    end

    end

    function [Label_Mask,Type_Mask,N_Obj] = sliderAnnotate(Img_Big_Seg_temp,Img_Big_Full)
        Label_Mask = zeros(size(Img_Big_Full));
        Type_Mask  = zeros(size(Img_Big_Full));
        %%% 1) Black 2) White 3) Red 4) Green
        Mycolormap  = [[0 0 0];[1 1 1];[1 0 0];[0 1 0]];
        fig = uifigure;
        fig.Position = [200 0 1200 800];
        ax     = axes(fig);
        ax.Position = [0.120 0.1100 0.7750 0.8150];
        Printing_Mask = double(Img_Big_Seg_temp);
        [Img_Big_Seg, N_objects] = bwlabel(Img_Big_Seg_temp,8);
        Obj_ptr = 1;
        Del_flag = 0;                          
        N_Obj = 1;
        Current_Mask = (Img_Big_Seg  == Obj_ptr);
        Printing_Mask(Current_Mask) = 3;
        Img_Ptr = imagesc(ax,Printing_Mask);
        colorarr = colorcube(40);
        str_list = ["Background";"Unlabeled";"Deleted";"Current";"Label 1"];
        N_labels = length(str_list);
        Mycolormap_temp  = [Mycolormap;colorarr(1,:)];
        colormap(ax,Mycolormap_temp);
        
        colorbar(ax,'Ticks',1:N_labels,'Color','Black','FontSize',17,'TickLabels',str_list);
        clim(ax,[0,N_labels]);
        set(ax,'XColor', 'none','YColor','none')

        uicontrol('Parent',fig,'Style','text','String',"Label",'Position',[30 680 100 20]);
        listbox_ptr_label = uicontrol('Parent',fig,'Style','listbox','String',["Delete","1"],'Position',[30 380 100 300]);
        uicontrol('Parent',fig,'Style','text','String',"Type",'Position',[30 300 100 20]);
        listbox_ptr_type = uicontrol('Parent',fig,'Style','listbox','String',["Soma","Axon","Background","Unkown?"],'Position',[30 200 100 100]);
        btn_back = uicontrol('Parent',fig,'Style', 'pushbutton','String',"Back"  ,'Position', [350 10 100 35],'Callback',@FnBackBtn,'Enable','off');
        uicontrol('Parent',fig,'Style', 'pushbutton','String',"Next"  ,'Position', [450 10 100 35],'Callback',@FnNextBtn);
        uicontrol('Parent',fig,'Style', 'pushbutton','String',"Add Label"  ,'Position', [30 330 100 35],'Callback',@AddLabel);
        listbox_ptr_label.Value = 2;
        uiwait(fig)

        function AddLabel(~,~)
            lst_str = string(listbox_ptr_label.String);
            N = length(lst_str);
            if (N<20)
                listbox_ptr_label.Value = 1:N;
                listbox_ptr_label.String = ["Delete",string(1:N)];
                str_list = [str_list;strcat("L   abel ",num2str(N))];
                N_labels = length(str_list);
                Mycolormap_temp  = [Mycolormap;colorarr(1:N,:)];
                colormap(ax,Mycolormap_temp);
                colorbar(ax,'Ticks',1:N_labels,'Color','Black','FontSize',17,'TickLabels',str_list);
                clim(ax,[0,N_labels]);
            end
        end

        function FnBackBtn(~,~)

            Current_Mask = (Img_Big_Seg == Obj_ptr);
            Previous_Mask = (Img_Big_Seg == Obj_ptr-1);
            Printing_Mask(Current_Mask) = 1;   %%% White
            Printing_Mask(Previous_Mask) = 3;   %%% Green
            Obj_ptr = Obj_ptr - 1;
            if  (Obj_ptr == 1)
                btn_back.Enable = 'off';
            end
            set(Img_Ptr,'CData', Printing_Mask);
        end

        function FnNextBtn(~,~)
            if (Obj_ptr == 1)
                btn_back.Enable  = 'on';
            end
            Current_Label = listbox_ptr_label.Value-1;
            Current_Type = listbox_ptr_type.Value;
            Current_Mask = (Img_Big_Seg == Obj_ptr);
            Label_Mask(Current_Mask) = Current_Label;
            Type_Mask(Current_Mask)  = Current_Type;

            if (Current_Label == 0)     %%% Printing what's done
                Printing_Mask(Current_Mask) = 2;   %%% Red
            else
                Printing_Mask(Current_Mask) = Current_Label+3;   %%% Blue
            end

            %%%% Step up pointer
            Obj_ptr = Obj_ptr + 1;
            if (Obj_ptr  <=  N_objects)     %%% Printing Current new
                Current_Mask  = (Img_Big_Seg == Obj_ptr);
                Printing_Mask(Current_Mask) = 3;      %%% Green
            end
             
            set(Img_Ptr,'CData', Printing_Mask);
            if (Obj_ptr > N_objects)
                FinishAnno();
            end
        end

        function FinishAnno(~,~)
            N_Obj = listbox_ptr_label.Value-1;
            delete(fig)
        end

    end
end
