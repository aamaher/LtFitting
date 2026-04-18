function Return_vector = Run_Sim_v5(path_main,folder_name_main,N_measure,Mode,Auto_threshold,param,param_plot)
    %% Salem - To do list:
    % 1- Add the autothreshold
    % 2- Test one frame code
    % 3- Add slidebar for min peak intensity
    % 4- Remove the flying checkbox, make objects filter always on
    % 5- Integrate registeration code here as well (Analysis_type = 3).
    %% Different modes available: (Masks, Labels, Results will be overwritten)
    % 1) Get Thresholded Mask
    % 2) Lifetime Analysis
    % 3) Annotate/Re-Annotate
    % 4) Plot results
    %% Setting folders and some parameters
    Return_vector{1} = struct("status",1);    
    warning('off','All')
    path_main = char(path_main);
    if (param.Analysis_type == 3 && param.ControlsCalib == 0)
        param.cell_level_analysis = 1;
        param.pixelwise = 0;
    else
        param.cell_level_analysis = 0;
    end

    if (param.ControlsCalib == 1 && param.Analysis_type == 3)
         param.pixelwise = 0;
    end
    if (path_main(end) == '/')
        path_main = path_main(1:end-1);
    end
    folder_name_main = char(folder_name_main);
    if (folder_name_main(end) == '/')
        folder_name_main = folder_name_main(1:end-1);
    end
    if (~isfolder('Results'))
        mkdir('Results')
    end
    if (~isfolder("Results/Masks"))
        mkdir('Results/Masks')
    end
    if (~isfolder("Results/AmpMap"))
        mkdir('Results/AmpMap')
    end

    %% Finding File names in the folder
    if (param.Analysis_type == 1)
        filename{1} = get_files(path_main,folder_name_main);
        if (N_measure == 0)
            N_measure  = 1:size(filename{1},2);            
        end
        sweep_count = length(N_measure);
        N_folders = 1;
        folder_name{1} = folder_name_main;
        aligned_files = filename{1};
        all_bases = filename{1};
        N_all_unique = 1;
    elseif (param.Analysis_type == 2 || param.Analysis_type == 3)
        path_main = fullfile(path_main,folder_name_main);
        files_temp = dir(path_main);
        dirFlags = [files_temp.isdir];
        subFolders = files_temp(dirFlags); % A structure with extra info.
        mainFiles = files_temp(~dirFlags); % A structure with extra info.        
        folder_name = {subFolders(3:end).name}; % Start at 3 to skip . and ..
        if (param.Stain == 1)
            mainFiles_names = {mainFiles.name}; % Start at 3 to skip . and ..
            mainFiles_names(endsWith(mainFiles_names, '.set')) = [];
            mainFiles_names(endsWith(mainFiles_names, '.txt')) = [];
        
            tok = cellfun(@(x) regexp(x, '_(\d+)_Filter\d+_\d{2}\.[^.]+$', 'tokens', 'once'), ...
                mainFiles_names, 'UniformOutput', false);
            stainfiles_idx = cellfun(@(x) str2double(x{1}), tok);
            if (isempty(stainfiles_idx))
                Check_stain_exist = 0;
            else
                Check_stain_exist = 1;
            end
        end
        folder_name(strcmp(folder_name, 'Properties')) = [];
        N_folders = length(folder_name);
        if isempty(param.PriorMask)
            param.PriorMask = zeros(1,N_folders);
            fprintf("\nparam.PriorMask was empty, filled by zeros\n");
        end
        if (isempty(param.FoldersToRun))
            param.FoldersToRun = 1:N_folders;
            fprintf("\nparam.FoldersToRun was empty (Analysis_type = 2), filled by zeros\n");
        end
        for idx_folder = 1:N_folders
            filename{idx_folder} = get_files(path_main,folder_name{idx_folder});
            size_arr(idx_folder) = size(filename{idx_folder},2);
        end
        for idx_folder = 1:N_folders
            str_arr_temp = string(filename{idx_folder}');
            %%% Removes Measurement Number "_01"
            str_arr_edit_temp = extractBefore(str_arr_temp, strlength(str_arr_temp) - 2);
            %%% Removes Filter/Day/WVL tags
            str_arr_edit_temp = regexprep(str_arr_edit_temp, '(?i)(Filter|Day|WVL|D)\d+', '');
            str_arr_idx_only{idx_folder} = double(regexprep(str_arr_edit_temp, '.*(?i)idx(\d+)_?.*', '$1'));
            %%% Removes underscores
            str_normalized{idx_folder} = regexprep(str_arr_edit_temp, '_', '');
            %%% Find unique base names ===
            str_unique{idx_folder} = unique(str_normalized{idx_folder});
        end

        
        % Find unique base names across folders ===
        all_bases = unique(vertcat(str_unique{:}));
        if (N_measure == 0)
            sweep_count = size(all_bases,1);
            N_measure  = 1:sweep_count;
        else
            sweep_count = length(N_measure);
        end

        N_all_unique = numel(all_bases);
        % Preallocate aligned filenames ===
        aligned_files = strings(N_folders, numel(all_bases));
        aligned_files(:) = missing;
        
        % === STEP 4: Fill in corresponding filenames ===
        for i = 1:N_folders
            [tf, loc] = ismember(str_normalized{i}, all_bases);
            aligned_files(i, loc(tf)) = filename{i}(tf); % keep original name with filter/day info
        end
        for i = 1:N_folders
            filename_new{i} = aligned_files(i, :);
        end
        filename = filename_new;
        % === STEP 5: Replace missing entries with NaN (optional) ===
        % aligned_files = replace(aligned_files, string(missing), "NaN");
        

       
        if isempty(param.PriorVector)
            param.PriorVector = zeros(1,N_folders);
        end
        if ismember(param.Analysis_type,[1,2])
            Return_vector = cell(N_folders,1);
            for idxx = 1:N_folders
                Return_vector{idxx} = struct("status",1);
            end
        end
        if (param.Analysis_type == 1)
            if (N_measure == 0)
                sweep_count = unique(size_arr);
                N_measure  = 1:sweep_count;
            else
                sweep_count = length(N_measure);
            end
        end
    end

    % === STEP 6 (optional): Make table for easy inspection ===
    if (param.DEBUGGING == 1)
        index_row = string(1:numel(all_bases));  % make numeric index labels
        T = array2table(aligned_files, ...
        'RowNames', folder_name, ...
        'VariableNames', cellstr(all_bases));

        T_index = array2table(index_row, 'VariableNames', cellstr(all_bases));
        T_index.Properties.RowNames = {'Index'};
        
        % Combine the index row with your main table vertically
        T_full = [T_index; T];
        
        disp(T_full);
    end

    
    if (param.Stain)
        if (Check_stain_exist)
            param.stain_exists_flag = ismember(str_arr_idx_only{1}, stainfiles_idx);
            param.stain_matching_files = arrayfun(@(k) mainFiles_names(stainfiles_idx == k), str_arr_idx_only{1},'UniformOutput', false);
        end
    end
    %% Mode 1 - Thresholding
    if (Mode == 1 || Mode == 0)
        arr_to_run = 1;
        if ~isempty(param.FoldersToRun)
            arr_to_run = param.FoldersToRun;
        end
        
        for j = 1:sweep_count
            idx = N_measure(j);
            if (param.Stain == 1)
                if (param.stain_exists_flag(idx))
                    param.current_stain_exists_flag = 1;
                    param.current_stain_filenames = param.stain_matching_files{idx};
                else
                    param.current_stain_exists_flag = 0;
                end
            end

            N_folders_to_Run = length(arr_to_run);
            path_files_temp = cell(N_folders_to_Run,1);
            Mask_filepath_arr = cell(N_folders_to_Run,1);
            selected_filenames = cell(N_folders_to_Run,1);
            for idx_folder_current = 1:N_folders_to_Run
                idx_folder = arr_to_run(idx_folder_current);
                current_filename = filename{idx_folder}(idx);
                if ismissing(current_filename)
                    filename_exist_arr(idx_folder_current) = false;
                    continue;
                end
                path_files_temp{idx_folder_current} = fullfile(path_main,folder_name{idx_folder},current_filename);
                selected_filenames{idx_folder_current} = current_filename;
                param.folder_idx_current = idx_folder;
                Mask_filepath_arr{idx_folder_current} = Get_savefile_name(param,path_files_temp{idx_folder_current},"Mask");
                filename_exist_arr(idx_folder_current) = true;
            end
            if (~param.Rerun_Mask)
                good_arr = arr_to_run(filename_exist_arr);
                N_good_folders = length(good_arr);
                for idx_folder_current = 1:N_good_folders
                    idx_good = good_arr(idx_folder_current);
                    if exist(Mask_filepath_arr{idx_good},'file')
                        flag_run = 0;
                    else
                        flag_run = 1;
                        break;
                    end
                end
            else
                flag_run = 1;
            end
            if (flag_run)
                fprintf("Thresholding: %d/%d (%d/%d)",j,sweep_count,idx,N_all_unique);
                for idx_folder_current = 1:N_folders_to_Run
                    fprintf(", %s",string(selected_filenames{idx_folder_current}));
                end
                fprintf("\n")
                if (param.Analysis_type == 1)
                    [segment_Obj,Selected_ROI] = ThresholdingSlider(param,path_files_temp{1},Auto_threshold,folder_name,selected_filenames);
                    Skip = 0;
                    arr_folders = 1;
                else
                    N_folders_to_Run = sum(filename_exist_arr);
                    arr_folders = find(filename_exist_arr == 1);
                    [segment_Obj,Skip] = ThresholdingSlider_Multi(param,path_files_temp,Auto_threshold,folder_name,selected_filenames,filename_exist_arr);
                    Selected_ROI = [];
                end
                Skip_once = 0;
                for idx_folder_current = 1:N_folders_to_Run 
                    idx_folder = arr_folders(idx_folder_current);
                    if Skip == 0
                        if (sum(sum(segment_Obj{idx_folder_current})) ~= 0)
                            Mask = (segment_Obj{idx_folder_current} ~= 0);
                        else
                            Mask = [];
                            Skip_once = 1;
                        end
                            
                    else
                        Mask = [];
                    end
                    save(Mask_filepath_arr{idx_folder},'Mask',"Selected_ROI",'Skip','Skip_once')
                end
            else
                if (param.FixOldMasks)
                    fprintf("Thresholding: %d/%d (%d/%d), Fixing old existing",j,sweep_count,idx,N_all_unique);
                    for idx_folder_current = 1:N_folders_to_Run
                        if (filename_exist_arr(idx_folder_current) == 0)
                            continue;
                        end
                        load(Mask_filepath_arr{idx_folder_current},'Mask',"Selected_ROI",'Skip')
                        Mask_OG{idx_folder_current} = Mask;
                            if isempty(Mask)
                                Skip_OG{idx_folder_current} = 1;
                            else
                                Skip_OG{idx_folder_current} = Skip; 
                            end
                    end                    
                    FixFilter(param,path_files_temp,filename_exist_arr,Mask_OG,Skip_OG);
                else
                    fprintf("Thresholding: %d/%d (%d/%d), Skipped since it exists",j,sweep_count,idx,N_all_unique);
                end
                for idx_folder_current = 1:N_folders_to_Run
                    fprintf(", %s",string(selected_filenames{idx_folder_current}));
                end
                fprintf("\n")
            end
        end
    end
    
    %% Mode 2 - Sweep across the Field of views and Perform the lifetime fitting for each one.
    if (Mode == 2 || Mode == 0 && Show_results_only == 0)
        if (param.parallel == 1 && param.DEBUGGING == 0)
            pool = gcp;             % Get the current parallel pool
            if ~isempty(pool)       % Parallel pool is active
                % fprintf('A parallel pool is active with %d workers.\n', pool.NumWorkers);
            else
                fprintf('No parallel pool is active.\n');
                parpool('threads')
            end
        end
        counter = 0;
        arr_to_run = 1;
        if ~isempty(param.FoldersToRun)
            arr_to_run = param.FoldersToRun;
        end
        Returned_data = cell(sweep_count,1);
        for j = 1:sweep_count
            idx = N_measure(j);
            PriorData = cell(N_folders,1);

            if (param.Stain == 1)
                if (param.stain_exists_flag(idx))
                    param.current_stain_exists_flag = 1;
                    param.current_stain_filenames = param.stain_matching_files{idx};
                else
                    param.current_stain_exists_flag = 0;
                end
            end

            if (param.Analysis_type ==3)
                path_list = cell(1, numel(param.FoldersToRun));
                skip_any  = false;
                for ii = 1:numel(param.FoldersToRun)
                    f = param.FoldersToRun(ii);
                    current_filename = filename{f};
                    if ismissing(current_filename(idx))
                        skip_any = true;
                    else
                        p = fullfile(path_main, folder_name{f}, current_filename{idx});
                        path_list{ii} = p;
                        Mask_filepath = Get_savefile_name(param, p, "Mask");
                        if isfile(Mask_filepath)
                            M = load(Mask_filepath);
                            Mask_arr{ii} = M.Mask;
                            if (isfield(M,'Skip') && M.Skip == 1) || (isfield(M,'Skip_once') && M.Skip_once == 1)
                                skip_any = true;
                            end
                        end
                    end
                end
        
                if skip_any
                    str_filename = current_filename(idx);
                    if ismissing(str_filename)
                        str_filename = "Missing";
                    end
                    fprintf('\n%d/%d, %d- %s. Skipped.\n', j, sweep_count, N_measure(j),str_filename);
                    continue;
                end
                if (~param.Rerun)
                    %%% Find me
                    filename_results_current = Get_savefile_name(param,path_list{1},"Results");
                    if exist(filename_results_current,'file')
                        continue;
                    end
                end
        
                [~, Returned_data{j}, status] = SalemPixelFitting22(path_list, Mask_arr, param);
                fprintf('%d/%d, %d- %s. Status: %d.\n', j, sweep_count, N_measure(j), current_filename{idx}, status);

            else
                Skip_all = 0;
                if (param.Mask_together == 1)
                    Mask_all = [];
                    for idx_folder = arr_to_run  %%% Salem: Take care of this
                        current_filename = filename{idx_folder};
                        if ismissing(current_filename(idx))
                            continue;
                        end
                        path = fullfile(path_main,folder_name{idx_folder},current_filename{idx});
                        Mask_filepath = Get_savefile_name(param,path,"Mask");
                        M = load(Mask_filepath);
                        if isempty(Mask_all)
                            Mask_all = M.Mask;
                        else
                            Mask_all = Mask_all.*M.Mask;
                        end
                    end
                    for idx_folder = arr_to_run  %%% Salem: Take care of this
                        Mask_arr{idx_folder} = Mask_all;
                    end
                    if isempty(Mask_all)
                        Skip_all = 1;
                    end
                end
                for idx_folder = arr_to_run  %%% Salem: Take care of this
                    param.PriorData = [];   %%% Reset for every folder
                    if (N_folders && ~isempty(param.order_map))
                        if length(param.order_map) == N_folders
                            param.order = param.order_map(idx_folder);
                        end
                    end
                    param.folder_idx_current = idx_folder;
                    counter = counter + 1;
                    current_filename = filename{idx_folder};
                    if ismissing(current_filename(idx))
                        continue;
                    end
                    path = fullfile(path_main,folder_name{idx_folder},current_filename{idx});
                    if (param.Mask_together == 1)
                        Mask = Mask_all;
                    else
                        Mask_filepath = Get_savefile_name(param,path,"Mask");
                        M = load(Mask_filepath);
                        Mask = M.Mask;
                    end
                    Selected_ROI = M.Selected_ROI;
                    if (isfield(M,'Skip') && M.Skip == 1) || (isfield(M,'Skip_once') && M.Skip_once == 1) || Skip_all
                        Skip = true;
                    else
                        Skip = false;
                    end
   
                    if (Skip == 1)
                        fprintf('\n%d/%d, %d-%s. Skipped \n',idx_folder+N_folders*(j-1),sweep_count*N_folders,N_measure(j),current_filename{N_measure(j)});
                        continue;
                    end
                    fprintf('%d/%d, %d-%s. ',idx_folder+N_folders*(j-1),sweep_count*N_folders,N_measure(j),current_filename{N_measure(j)});
                    param.Mask = Mask;
                    param.ROI = Selected_ROI;
                    if (param.Analysis_type == 2)
                        if (param.PriorVector(idx_folder) > 0)
                            prior_filename = filename{param.PriorVector(idx_folder)};
                            path_prior = fullfile(path_main,folder_name{param.PriorVector(idx_folder)},prior_filename{idx});
                            param.folder_idx_current = param.PriorVector(idx_folder);
                            Results_filepath = Get_savefile_name(param,path_prior,"Results1D");
                            param.folder_idx_current = idx_folder;
                            if(~isfile(Results_filepath))
                                disp("Can't find Fitting Results")
                                return
                            else
                                load(Results_filepath,'tau_map')
                            end
                            param.PriorData = tau_map;
                            if (length(param.order_map)>1)
                                param.order = param.order_map(idx_folder);
                            end
                        end
                    end
                    if (~param.Rerun)
                        filename_results_current = Get_savefile_name(param,path,"Results");
                        if exist(filename_results_current,'file')
                            fprintf(', Skipped since results found\n');        
                            continue;
                        end
                    end
                    
                    [PriorData{idx_folder},Returned_data, status] = SalemPixelFitting22(path,Mask_arr, param);
                    fprintf(', Status: %d.\n',status);        
                end
            end
        end
    Return_vector = Returned_data;
    end
    
    %% Mode 3 - Segmenting and labeling Objects
    if (Mode == 3)
        for j= 1:sweep_count
            for idx_folder = 1:N_folders
                idx = N_measure(j);
                param.folder_idx_current = idx_folder;
                keyboard;         %%% Salem: Make sure this is correct here
                path_current = fullfile(path_main,filename{idx});
                Mask_filepath = Get_savefile_name(param,path_current,"Mask");
                Label_filepath = Get_savefile_name(param,path_current,"Label");
                load(Mask_filepath,'Mask',"Selected_ROI")
                [Label_Mask,Type_Mask,N_Obj] = Obj_segmentation(path,param,param_plot,Mask);
                save(Label_filepath,'Label_Mask','Type_Mask','N_Obj','-mat');
            end
        end
    end
    
    %% Mode 4 - Plot the Results
    if (Mode == 4 || Mode == 0)
        arr_to_run = 1;
        for j = 1:sweep_count
            if ~isempty(param.FoldersToRun)
                arr_to_run = param.FoldersToRun;
            end
            if (param.Analysis_type == 3)
                arr_to_run = 1;  %%% We plot only results from one file named after the first folder
            end
            for idx_folder = arr_to_run
                if (length(param.order_map)>1)
                    param.order = param.order_map(idx_folder);
                end
                current_filename = filename{idx_folder};
                param.folder_idx_current = idx_folder;
                if (param.Analysis_type == 2)
                    if (param.PriorVector(idx_folder) > 0)
                        param.PriorCell = 1;
                    else
                        param.PriorCell = [];
                    end
                end
                idx = N_measure(j);
                path = fullfile(path_main,folder_name{idx_folder},current_filename(idx)); 
                if ismissing(current_filename(idx))
                    fprintf("\n%d - Missing",idx);
                    continue;
                else
                    fprintf("\n%d - %s",idx,current_filename{idx});
                end
                
                if (param.Stain == 1)
                    if (param.stain_exists_flag(idx))
                        param.current_stain_exists_flag = 1;
                        param.current_stain_filenames = param.stain_matching_files{idx};
                    else
                        param.current_stain_exists_flag = 0;
                    end
                end

                Returned_data = plotObjLabels11_Labels(param,param_plot,path,Return_vector{idx_folder});
                if (param.Analysis_type == 3)
                    Return_vector{j} = Returned_data;
                else
                    Return_vector{idx_folder} = Returned_data;
                end
            end
        end
    end
end

%% File loading function
function filename = get_files(path_main,folder_name)
    path_current = fullfile(path_main,folder_name);
    path_ls = dir(path_current);
    j = 1;
    filename = {};
    for i = 3:size(path_ls,1)
        if (path_ls(i).isdir == 1)
            if (path_ls(i).name == "RecSettings.txt" || path_ls(i).name == "Properties")
                continue;
            end    
            filename{j} = path_ls(i).name;
            j = j+1;
        end
    end   
end

%% Thresholding Function 
function [Mask,Selected_ROI] = ThresholdingSlider(param,path,Auto_threshold,folder_name,selected_filenames)
    path     = char(path);
    Selected_ROI = [];
    N_folder = size(path,1);
    I = cell(N_folder,1);
    for idx_folder = 1:N_folder
        files = dir([path(idx_folder,:), '/*.tif']);
        parametersFile = strcat(path(idx_folder,:),'/RecSettings.txt'); % contains experimental parameters
        [~,delta_t,tmin,tmax,~,BINNING,~,~] = GetExperimentalParameters(parametersFile);    % get experimental parameters
        t_sig = (tmin:delta_t:tmax)';
        n = length(t_sig);
        image_filenames = sort(split(strip(ls([path(idx_folder,:), '/*.tif']))));   % image files (sorted)
        HorizontalPixels = 1:1936/BINNING;     % based on binning factor [pixels]
        VerticalPixels = 1:1216/BINNING;         % based on binning factor [pixels]
        if isempty(param.idx_peak_thresh)
            I_temp = zeros(range(VerticalPixels)+1,range(HorizontalPixels)+1,n);       % matrix which array C will be converted to
            for i = 1:n     %% extract pixel intensities and time delays for each frame
                image_filename_temp = erase(string(image_filenames{i}) , "'" );
                Image_handle = Tiff(image_filename_temp,'r');  % read given frame
                I_temp(:,:,i) = read(Image_handle); % extract intensity data for given frame
            end
            I_sum = squeeze(sum(I_temp,[1,2]));
            [~,idx_peak] = max(I_sum);
            I{idx_folder} = I_temp(:,:,idx_peak);
        else
            image_filename_temp = erase(string(image_filenames{param.idx_peak_thresh}) , "'" );
            Image_handle = Tiff(image_filename_temp,'r');  % read given frame
            I{idx_folder} = double(read(Image_handle));
        end
        if (Auto_threshold ~= 0)
            Threshold = Auto_threshold;
            Img = imsharpen(I{idx_folder},'Radius',30,'Amount',2);
            peak = max(Img(:));
            Mask_temp = (Img > Threshold*peak);
            [Mask_temp, ~] = bwlabel(Mask_temp ,8);
            Mask_temp(I{idx_folder}<400) = 0;
            Mask{idx_folder}  = bwareaopen(Mask_temp ,50,8);
        end
    end
    if (length(param.PriorMask) == N_folder)
        Mask_common = Mask{1};
        for idx_folder = 2:N_folder
            Mask_common = Mask_common.*Mask{idx_folder};
        end
        for idx_folder = 1:N_folder
            Mask{idx_folder} = Mask_common;
        end
    else
        PriorMask = [];
    end
    if(Auto_threshold == 0)
        [Mask,Selected_ROI] = sliderImageDemo(I,folder_name,selected_filenames,PriorMask);
    end
    
end

%% Slider Main
function [Mask_calc,ROI] = sliderImageDemo(Peak_Imgs,folder_name,selected_filenames,PriorMask)
    N_folder = size(Peak_Imgs,2);
    flag_ROI = 0;
    flag_sharpness = 1;
    flag_filter = 1;
    ROI = [];
    fig = uifigure;
    fig.Name = "Setting Mask";
    if (N_folder == 1)
        fig.Position = [50, 50, 1400, 500];
    else
        fig.Position = [50, 50, 1400, 750];
    end
    fig.AutoResizeChildren = 'off';
    X = size(Peak_Imgs,1)/2;
    Y = size(Peak_Imgs,1)/2;
    L = 100;
    W = 100;
    h = -30;
    w = 150;
    w1 = 150;
    Thresh = cell(N_folder,1);
    N_obj = cell(N_folder,1);
    A = cell(N_folder,1);
    R = cell(N_folder,1);
    None_flag = cell(N_folder,1);
    Union_flag = cell(N_folder,1);
    Intersect_flag = cell(N_folder,1);
    idx_folder = 1;
    Thresh{idx_folder} = 0.05;
    N_obj{idx_folder} = 150;
    A{idx_folder} = 2;
    R{idx_folder} = 30;

    if (~isempty(PriorMask))
        for idx_folder  = 1:N_folder
            Thresh{idx_folder} = 0.05;
            N_obj{idx_folder} = 150;
            A{idx_folder} = 2;
            R{idx_folder} = 30;
            if (PriorMask(idx_folder) == 0)
                None_flag{idx_folder} = 1;
                Union_flag{idx_folder} = 0;
                Intersect_flag{idx_folder}  = 0;
            elseif (PriorMask(idx_folder) > 0)
                None_flag{idx_folder} = 0;
                Union_flag{idx_folder} = PriorMask(idx_folder);
                Intersect_flag{idx_folder}  = 0;
            elseif (PriorMask(idx_folder) < 0)
                None_flag{idx_folder} = 0;
                Union_flag{idx_folder} = 0;
                Intersect_flag{idx_folder} = abs(PriorMask(idx_folder));
            end
        end
    end

    Thresh_OG = Thresh;
    N_obj_OG = N_obj;
    A_OG = A;
    R_OG = R;
    if (N_folder == 1)
        frame = uipanel(fig,'Position',[w+900, 20, 300, 470]);
    else
         frame = uipanel(fig,'Position',[w+900, 10, 300, 700]);
    end
    frame.Title = 'Controls';
    %%% Frames
    frame_T = uipanel(frame,"Position",[25, 365, 155, 80],"Title",strcat("Threshold: ",num2str(Thresh{1})));
    frame_Sharpness = uipanel(frame,"Position",[25, 275, 250, 80],"Title","Sharpness");
    frame_ROI = uipanel(frame,"Position",[25, 105, 250, 160],"Title","Custom ROI");
    frame_N = uipanel(frame,"Position",[25, 40, 250, 60],"Title","Obj. Filer");

    %%% Threshold frame elements
    slider_T = uicontrol('Parent',frame_T,'Style', 'slider', 'Min', 0, 'Max', 0.8, 'Value', Thresh{idx_folder}, 'Position', [10, 20, 130, 20], 'Callback', @UpdateThresh);

    %%% Sharpness frame elements
   
    slider_A = uicontrol('Parent',frame_Sharpness,'Style', 'slider','Position', [55, 30, 100, 20], 'Min', 0, 'Max', 2, 'Value', A{idx_folder}, 'Callback', @UpdateThresh);
               uicontrol('Parent',frame_Sharpness,'Style', 'text',  'Position', [30, 30, 20, 20],'String', 'A');
    label_A =  uicontrol('Parent',frame_Sharpness,'Style', 'text',  'Position', [155, 30, 30, 20],'String', num2str(A{idx_folder}));
    slider_R = uicontrol('Parent',frame_Sharpness,'Style', 'slider','Position', [55, 10, 100, 20], 'Min', 0, 'Max', 30, 'Value', R{idx_folder},'Callback', @UpdateThresh);
               uicontrol('Parent',frame_Sharpness,'Style', 'text',  'Position', [30, 5, 20, 20],'String', 'R');
    label_R =  uicontrol('Parent',frame_Sharpness,'Style', 'text',  'Position', [155, 5, 30, 20],'String', num2str(R{idx_folder}));
    checkbox_Sharpness = uicheckbox('Parent',frame_Sharpness,'Position',[10, 35, 100, 22],'Text',' ','Value',1,'ValueChangedFcn',@toggleSharpness);
    
    % ROI
    checkbox = uicheckbox(frame_ROI, 'Position', [10, 110, 100, 22],'ValueChangedFcn',@toggleSliderVisibility);
    slider_L = uicontrol('Parent',frame_ROI,'Style', 'slider', 'Min', 0, 'Max', size(Peak_Imgs{1},2), 'Value', L, 'Position', [60, 85, 100, 20], 'Callback', @updateParams);
    slider_W = uicontrol('Parent',frame_ROI,'Style', 'slider', 'Min', 0, 'Max', size(Peak_Imgs{1},1), 'Value', W, 'Position', [60, 60, 100, 20], 'Callback', @updateParams);
    slider_Y = uicontrol('Parent',frame_ROI,'Style', 'slider', 'Min', 0, 'Max', size(Peak_Imgs{1},1), 'Value', Y, 'Position', [60, 35, 100, 20], 'Callback', @updateParams);
    slider_X = uicontrol('Parent',frame_ROI,'Style', 'slider', 'Min', 0, 'Max', size(Peak_Imgs{1},2), 'Value', X, 'Position', [60, 10, 100, 20], 'Callback', @updateParams);

    uicontrol('Parent',frame_ROI,'Style', 'text', 'Position', [10, 85, 50, 20], 'String', 'Length');
    uicontrol('Parent',frame_ROI,'Style', 'text', 'Position', [10, 60, 50, 20], 'String', 'Width');
    uicontrol('Parent',frame_ROI,'Style', 'text', 'Position', [10, 35, 50, 20], 'String', 'Y');
    uicontrol('Parent',frame_ROI,'Style', 'text', 'Position', [10, 10, 50, 20], 'String', 'X');

    label_L = uicontrol('Parent',frame_ROI,'Style', 'text', 'Position', [160, 85, 50, 20],'String', num2str(L));
    label_W = uicontrol('Parent',frame_ROI,'Style', 'text', 'Position', [160, 60, 50, 20],'String', num2str(W));
    label_Y = uicontrol('Parent',frame_ROI,'Style', 'text', 'Position', [160, 35, 50, 20],'String', num2str(Y));
    label_X = uicontrol('Parent',frame_ROI,'Style', 'text', 'Position', [160, 10, 50, 20],'String', num2str(X));


    %%% Objects filter
    
    slider_N = uicontrol('Parent',frame_N ,'Style', 'slider', 'Min', 0, 'Max', 200, 'Value', N_obj{idx_folder},'Position', [30, 10, 80, 20], 'Callback', @UpdateThresh);
    label_N = uicontrol('Parent',frame_N,'Style', 'text', 'Position', [110, 10, 30, 20], 'String', num2str(N_obj{idx_folder}));
    checkbox3 = uicheckbox(fig, 'Text', '  ','Position',[10, 10, 10, 22],'Value',1,'ValueChangedFcn',@toggleFilter);

    %%% Enable initially
    slider_L.Enable = 'off';
    slider_W.Enable = 'off';
    slider_X.Enable = 'off';
    slider_Y.Enable = 'off';
    slider_N.Enable = 'on';
    slider_A.Enable = 'on';
    slider_R.Enable = 'on';
    
    uicontrol('Parent',frame,'Style', 'pushbutton','String',"Start Analyzing",'Position', [50, 5, 100, 30],'Callback',@SetT)
    
    
    if (N_folder >1)
        bg_member = cell(N_folder,1);
        bg_mask =  cell(N_folder,1);
        bg = uibuttongroup(fig,"Position",[w1+925, h+500, 80, 30*N_folder],"Title","Edit Mask","SelectionChangedFcn",@SwitchMask);
        frame_register = uipanel(fig,"Position",[w1+1095, h+405, 80, 80],"Title","Register");
        uicontrol('Parent',frame_register,'Style', 'pushbutton','String',">",'Position', [50, 20, 20, 20],'Enable','off')
        uicontrol('Parent',frame_register,'Style', 'pushbutton','String',"<",'Position', [10, 20, 20, 20],'Enable','off')
        uicontrol('Parent',frame_register,'Style', 'pushbutton','String',"v",'Position', [30, 5, 20, 20],'Enable','off')
        uicontrol('Parent',frame_register,'Style', 'pushbutton','String',"^",'Position', [30, 35, 20, 20],'Enable','off')

        bg_process_relation = uibuttongroup(fig,"Position",[w1+1095, h+500, 80, 30*N_folder],"Title","Relation","SelectionChangedFcn",@ProcessChange);
        for idx_folder = 1:N_folder
            bg_mask{idx_folder} = uibuttongroup(fig,"Position",[w1+840+idx_folder*85, h+600, 80, 30*N_folder],"Title",strcat("Mask: ",num2str(idx_folder)),"SelectionChangedFcn",@updateImage);
            bg_process_relation_member{idx_folder} = uiradiobutton(bg_process_relation,"Text",strcat("Mask: ",num2str(idx_folder)),"Position",[10 50-(idx_folder-1)*20 100 15],'Value',1);
        end
        bg_process = uibuttongroup(fig,"Position",[w1+1010, h+500, 80, 30*N_folder],"Title","Processes","SelectionChangedFcn",@ProcessChange);
        bg_none = uiradiobutton(bg_process,"Text","None","Position",[10 45 100 22],'Value',1);
        bg_union = uiradiobutton(bg_process,"Text","Union","Position",[10 25 100 22],'Value',0);
        bg_intersect = uiradiobutton(bg_process,"Text","Intersect","Position",[10 5 100 22],'Value',0);
        bg_member_select = cell(N_folder,N_folder);
    end

    ax = cell(N_folder,1);
    bx = cell(N_folder,1);
    Mask_calc = cell(N_folder,1);
    imgHandle1 = cell(N_folder,1);
    rectangleHandle1 = cell(N_folder,1);
    rectangleHandle2 = cell(N_folder,1);
    for idx_folder = 1:N_folder
        ax{idx_folder} = axes('Parent', fig);
        bx{idx_folder} = axes('Parent', fig);
        if (N_folder == 1)
            ax{idx_folder}.Position = [0.05 0.25 0.3 0.5];
            bx{idx_folder}.Position = [0.4 0.25 0.3 0.5];
        else
            ax{idx_folder}.Position = [0.035+(idx_folder-1)*0.25 0.6 0.2 0.35];
            bx{idx_folder}.Position = [0.035+(idx_folder-1)*0.25 0.15 0.2 0.35];
            bg_member{idx_folder} = uiradiobutton(bg,"Text",strcat("Mask: ",num2str(idx_folder)),"Position",[10 45-(idx_folder-1)*20 100 22]);
            for idx_folder_select = 1:N_folder
                bg_member_select{idx_folder,idx_folder_select} = uiradiobutton(bg_mask{idx_folder},"Text",strcat("Mask: ",num2str(idx_folder_select)),"Position",[10 45-(idx_folder_select-1)*20 100 22]);
                if (idx_folder == idx_folder_select)
                    bg_member_select{idx_folder,idx_folder_select}.Value = 1;
                end
            end
        end
        Mask_calc_OG{idx_folder} = Thresh_Image(Thresh{idx_folder},idx_folder,Peak_Imgs{idx_folder});
    end
    
    for idx_folder = 1:N_folder
        if (N_folder > 1)
            if (PriorMask(idx_folder) == 0)
                if (idx_folder == 1)
                    bg_process_relation.SelectedObject = bg_process_relation_member{1};
                end
                Mask_calc{idx_folder} = Mask_calc_OG{idx_folder};
            else
                if (idx_folder == 1)
                    bg_process_relation.SelectedObject = bg_process_relation_member{abs(PriorMask(idx_folder))};
                end
                if (PriorMask(idx_folder) > 0)
                    if (idx_folder == 1)
                        bg_process.SelectedObject = bg_union;
                    end
                    Mask_calc{idx_folder} = Mask_calc_OG{idx_folder} | Mask_calc_OG{PriorMask(idx_folder)};
                elseif (PriorMask(idx_folder) < 0)
                    if (idx_folder == 1)
                        bg_process.SelectedObject = bg_intersect;
                    end
                    Mask_calc{idx_folder} = Mask_calc_OG{idx_folder} & Mask_calc_OG{abs(PriorMask(idx_folder))};
                end
            end
        else
            Mask_calc{idx_folder} = Mask_calc_OG{idx_folder};
        end
        if (flag_filter == 1)
            Mask_calc{idx_folder} = bwareaopen(Mask_calc{idx_folder} ,N_obj{idx_folder},8);
        end 
        Mask_temp = Mask_calc{idx_folder};
        Mask_temp(Peak_Imgs{idx_folder} < 400) = 0;
        Mask_calc{idx_folder} = Mask_temp;
        imgHandle1{idx_folder} = imagesc(bx{idx_folder},Mask_calc_OG{idx_folder});
        rectangleHandle1{idx_folder} = rectangle('Parent',ax{idx_folder},'Position', [X, Y, W, L], 'EdgeColor', 'b', 'FaceColor', 'none','Visible','off');
        rectangleHandle2{idx_folder} = rectangle('Parent',bx{idx_folder},'Position', [X, Y, W, L], 'EdgeColor', 'b', 'FaceColor', 'none','Visible','off');

        colormap(bx{idx_folder} , gray(2));
        max_arr = maxk(Peak_Imgs{idx_folder}(:),10);
        max_arr = mean(max_arr(2:end));
        N_pixles = sum(Mask_calc{idx_folder}(:));
        title(strcat("Mask: ",num2str(idx_folder),", Pixels: ",num2str(N_pixles),",Peak: ",num2str(max_arr)),'Parent',bx{idx_folder});
        Img_overlay = labeloverlay(imsharpen(imadjust(Peak_Imgs{idx_folder}./max(max(Peak_Imgs{idx_folder})))),Mask_calc{idx_folder},'Colormap',jet,'Transparency',0.1);
        imgHandle0{idx_folder} = imagesc(ax{idx_folder},Img_overlay);
        colormap(ax{idx_folder}, gray(256));
        title(strcat(string(folder_name{idx_folder}),': ',selected_filenames{idx_folder}),'Parent',ax{idx_folder},'Interpreter', 'none');
        if (N_folder == 1)
            cb = colorbar(bx{idx_folder});
            clim(bx{idx_folder},[0,0.7]);
            cb.Ticks = [1,2]; 
            cb.YTick = cb.Limits;
            cb.TickLabels = {'Off','On'};
        end
        ax{idx_folder}.UserData = linkprop([ax{idx_folder},bx{idx_folder}],{'x','y','xlim','ylim'});
        if (idx_folder > 1)
            ax{idx_folder}.UserData = linkprop([ax{idx_folder},ax{1}],{'x','y','xlim','ylim'});
            bx{idx_folder}.UserData = linkprop([bx{idx_folder},ax{1}],{'x','y','xlim','ylim'});
        end
    end
    if (N_folder >1)
        ax_temp = axes('Parent', fig);
        ax_temp.InnerPosition = [0.285,0.01,0.2,0.3];
        colormap(ax_temp,gray(2))
        cb = colorbar(ax_temp,"southoutside");
        cb.Ticks = [1,2]; 
        cb.YTick = cb.Limits;
        cb.TickLabels = {'Off','On'};
        cb.FontSize = 15 ;
        ax_temp.Visible = 'off';
    end
    uiwait(fig)

    
    function UpdateThresh(~,~)
        updateParams
        updateImage
    end
    
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
            for idx_folder = 1:N_folder
                rectangleHandle1{idx_folder}.Visible = 'on';
                rectangleHandle2{idx_folder}.Visible = 'on';
            end
        else
            slider_L.Enable  = 'off';
            slider_W.Enable  = 'off';
            slider_X.Enable  = 'off';
            slider_Y.Enable  = 'off';
            flag_ROI = 0;
            for idx_folder = 1:N_folder
                rectangleHandle1{idx_folder}.Visible = 'off';
                rectangleHandle2{idx_folder}.Visible = 'off';
            end
        end
    end
        
    function ProcessChange(~,~)
        if bg_process_relation.SelectedObject == bg_process_relation_member{1}
            V = 1;
        elseif bg_process_relation.SelectedObject == bg_process_relation_member{2}
            V = 2;
        elseif bg_process_relation.SelectedObject == bg_process_relation_member{3}
            V = 3;
        end
        for idx_folder = 1:N_folder
            if (bg_member{idx_folder}.Value ~= 1)
                continue;
            else
                Intersect_flag{idx_folder} = V*bg_intersect.Value;
                Union_flag{idx_folder} = V*bg_union.Value;
                None_flag{idx_folder} = V*bg_none.Value;
            end
        end
        updateImage;
    end

    function toggleSharpness(~,~)
            if checkbox_Sharpness.Value
                slider_A.Enable  = 'on';
                slider_R.Enable  = 'on';
                flag_sharpness = 1;
                A{idx_folder_current} = get(slider_A, 'Value');
                R{idx_folder_current} = get(slider_R, 'Value');
            else
                flag_sharpness = 0;
                slider_A.Enable  = 'off';
                slider_R.Enable  = 'off';
                set(imgHandle0{idx_folder}, 'CData', Peak_Imgs{idx_folder});
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
            set(imgHandle0{idx_folder}, 'CData', Peak_Imgs);
        end
        updateImage
    end

    function SwitchMask(~,~)
        for idx_folder_current = 1:N_folder
            if (bg_member{idx_folder_current}.Value ~= 1)
                continue;
            else
                set(frame_T, 'Title', strcat("Threshold: ",num2str(round(Thresh_OG{idx_folder_current},3))));
                set(label_A, 'String', round(A_OG{idx_folder_current},2));
                set(label_R, 'String', round(R_OG{idx_folder_current},2));
                set(label_N, 'String', round(N_obj_OG{idx_folder_current},2));

                set(slider_T, 'Value', round(Thresh_OG{idx_folder_current},3));
                set(slider_A, 'Value', round(A_OG{idx_folder_current},2));
                set(slider_R, 'Value', round(R_OG{idx_folder_current},2));
                set(slider_N, 'Value', round(N_obj_OG{idx_folder_current},2));
                
                if (None_flag{idx_folder_current} > 0)
                    bg_process.SelectedObject = bg_none;
                    bg_process_relation.SelectedObject = bg_process_relation_member{idx_folder_current};
                elseif (abs(Intersect_flag{idx_folder_current}) > 0)
                    bg_process.SelectedObject = bg_intersect;
                    V = abs(Intersect_flag{idx_folder_current});
                    bg_process_relation.SelectedObject = bg_process_relation_member{V};
                elseif (abs(Union_flag{idx_folder_current}) > 0)
                    bg_process.SelectedObject = bg_union;
                    V = Union_flag{idx_folder_current};
                    bg_process_relation.SelectedObject = bg_process_relation_member{V};
                end
                
            end
        end
    end
    
    function updateParams(~,~)
        %%% This updates the OG parameters. These parameters are only
        %%% updated if the OG mask is selected
        for idx_folder = 1:N_folder
            if (exist('bg_member','var'))
                V = bg_member{idx_folder}.Value;
            else
                V = 1;
            end
            if (V ~= 1)
                continue;
            else
                Thresh_OG{idx_folder} = get(slider_T, 'Value');    
                if (flag_filter == 1)
                    N_obj_OG{idx_folder} = round(get(slider_N, 'Value'));
                    if N_folder == 1    %%%% Salem: A quick fix, needs to be checked later
                        N_obj{idx_folder} = N_obj_OG{idx_folder};
                    end
                    set(label_N,'String',num2str(N_obj_OG{idx_folder}));
                end
                A_OG{idx_folder} = get(slider_A, 'Value');
                R_OG{idx_folder} = get(slider_R, 'Value');
            end
            Mask_calc_OG{idx_folder} = Thresh_Image(Thresh_OG{idx_folder},idx_folder,Peak_Imgs{idx_folder});
            set(imgHandle1{idx_folder}, 'CData', Mask_calc_OG{idx_folder});
            set(frame_T, 'Title', strcat("Threshold: ",num2str(round(Thresh_OG{idx_folder},3))));
            set(label_A, 'String', round(A_OG{idx_folder},2));
            set(label_R, 'String', round(R_OG{idx_folder},2));
        end
        
        if (flag_ROI == 1)
            L = round(get(slider_L, 'Value'));
            W = round(get(slider_W, 'Value'));
            X = round(get(slider_X, 'Value'));
            Y = round(get(slider_Y, 'Value'));
            for idx_folder = 1:N_folder
                rectangleHandle1{idx_folder}.Position = [X,Y,W,L];
                rectangleHandle2{idx_folder}.Position = [X,Y,W,L];
            end
        end
        set(label_L, 'String', round(L,2));
        set(label_W, 'String', round(W,2));
        set(label_X, 'String', round(X,2));
        set(label_Y, 'String', round(Y,2));
    end

    function updateImage(~,~)
        for idx_folder = 1:N_folder
            for idx_folder_ij = 1:N_folder
                if (exist('bg_member_select','var'))
                    V = bg_member_select{idx_folder,idx_folder_ij}.Value;
                else
                    V = 1;
                end
                if (V ~= 1)
                    continue;
                else
                    if (None_flag{idx_folder} > 0)
                        Mask_calc{idx_folder} = Mask_calc_OG{idx_folder_ij};
                    elseif(Intersect_flag{idx_folder} > 0)
                        Mask_calc{idx_folder} = Mask_calc_OG{idx_folder_ij} & Mask_calc_OG{Intersect_flag{idx_folder}};
                    elseif(Union_flag{idx_folder} > 0)
                        Mask_calc{idx_folder} = Mask_calc_OG{idx_folder_ij} | Mask_calc_OG{Union_flag{idx_folder}};
                    end
                    if (flag_filter == 1)
                        Mask_calc{idx_folder} = bwareaopen(Mask_calc_OG{idx_folder} ,N_obj{idx_folder},8);
                    end 
                    Mask_temp = Mask_calc{idx_folder};
                    Mask_temp(Peak_Imgs{idx_folder} < 400) = 0;
                    Mask_calc{idx_folder} = Mask_temp;
                    Img_overlay = labeloverlay(imsharpen(imadjust(Peak_Imgs{idx_folder}./max(max(Peak_Imgs{idx_folder})))),Mask_calc{idx_folder},'Colormap',jet,'Transparency',0.1);
                    set(imgHandle0{idx_folder}, 'CData', Img_overlay);
                    N_pixles = sum(Mask_calc{idx_folder}(:));
                    title(strcat("Mask: ",num2str(idx_folder),", Pixels: ",num2str(N_pixles),",Peak: ",num2str(max_arr)),'Parent',bx{idx_folder});
                end
            end
        end
    end

    function SetT(~,~)
        delete(fig)
    end

    function Mask = Thresh_Image(T_input,idx_folder_current,Image_input)    
        if (flag_sharpness == 1)
            A{idx_folder_current} = get(slider_A, 'Value');
            R{idx_folder_current} = get(slider_R, 'Value');
            Image = imsharpen(Image_input,'Radius',R{idx_folder_current},'Amount',A{idx_folder_current});
        end
        peak = max(Image(:));
        Mask = (Image > T_input*peak);
        Mask(Peak_Imgs{idx_folder_current} < 400) = 0;  %%% Delete Pixels with peak value less than 100 to avoid fitting errors
        if (flag_filter == 1)
            Mask  = bwareaopen(Mask ,N_obj{idx_folder_current},8);
        end 
    end
    
end

%% Segmentation Function
function [Label_Mask,Type_Mask,N_Obj] = Obj_segmentation(path,param,param_plot,Mask)
    
    tau_map = [];
    chi_map = [];
    filename_results_1D = Get_savefile_name(param,path,"Results1D");
    if (isfile(filename_results_1D))
        load(filename_results_1D,'tau_map','chi_map')
    else
        disp("No Mono-exponential (discrete) Lifetime data")
    end

    path_current   = char(path);
    if (path_current(end) == '/')
        path_current = path_current(1:end-1);
    end
    files = dir([path_current, '/*.tif']);
    image_filename =  fullfile(path_current,'/',files(param.idx_peak_thresh).name);
    Image_handle = Tiff(image_filename,'r');  % read given frame
    I_Peak = double(read(Image_handle)); % extract intensity data for given frame
    
    %%%  Using Watershed Segmentation Method
    bw  = Mask;
    D = -bwdist(~bw,'quasi-euclidean');
    Ld = watershed(D);
    bw2 = bw;
    bw2(Ld == 0) = 0;
    mask_em = imextendedmin(D,2);
    D2 = imimposemin(D,mask_em);
    Ld2 = watershed(D2);
    bw3 = bw;
    bw3(Ld2 == 0) = 0;
    Mask_Seg = ~bwareaopen(~bw3, 40);
    
%%  Annotate
    [Label_Mask,Type_Mask,N_Obj] = sliderAnnotate(Mask_Seg,I_Peak,tau_map,chi_map,param_plot);

%% Plot Intensity image on top of the intensity image
    function [Label_Mask,Type_Mask,N_Obj] = sliderAnnotate(Mask,I_peak,tau_map,chi_map,param_plot)
        
        fig = uifigure;
        ax     = axes(fig);
        ax2     = axes(fig);
        ax2_2     = axes(fig);
        ax2_temp     = axes(fig);
        ax3     = axes(fig);
        ax3_temp     = axes(fig);

        fig.Position = [150 0 1400 800];
        ax.Position = [0.1 0.1500 0.5 0.7];
        ax3.Position = [0.7 0.100 0.2 0.3];
        ax2.Position = [0.7 0.550 0.2 0.3];
        ax3_temp.Position = [0.55 0.1 0.4 0.3];
        ax2_temp.Position = [0.55 0.55 0.4 0.3];
        
        Label_Mask = zeros(size(I_peak));
        Type_Mask  = zeros(size(I_peak));
        
        Mycolormap  = [[0 0 0];[1 1 1];[1 0 0];[0 1 0]];    %%% 1) Black 2) White 3) Red 4) Green
        %%% Peak Intensity figure
        xmin = 200;
        xmax = 700;
        imagesc(ax2,I_peak/max(I_peak(:)));
        clim(ax2,[0,0.6]);
        xlim(ax2,[xmin,xmax]);
        colormap(ax2,'gray')
        if (~isempty(tau_map))
            %%% Lifetime map figure
            imagesc(ax2_2,tau_map,'AlphaData', 1*Mask);
            title(ax2,"Lifetime Map",'FontSize',15)
            ax2_2.UserData = linkprop([ax2,ax2_2],{'x','y','Position','xlim'});
            colormap(ax2_2,'turbo')
            xlim(ax2_2,[xmin,xmax]);
            clim(ax2_2,[param_plot.clim_tau_min,param_plot.clim_tau_max])    

            %%% Colorbar
            colormap(ax2_temp,'turbo')
            cb = colorbar(ax2_temp);
            cb.Title.String = '\tau (ns)';
            cb.FontSize = 15 ; 
            clim(ax2_temp,[param_plot.clim_tau_min,param_plot.clim_tau_max])

            %%% Chi map figure
            imagesc(ax3,chi_map);
            title(ax3,"\chi^{2} Map",'FontSize',15)            
            colormap(ax3,'turbo')
            xlim(ax3,[xmin,xmax]);

            hbins_chi = param_plot.hbins_chi;
            chi_avg = mean(chi_map(chi_map>0));
            chi_map (chi_map == 0) = NaN;
            [V_chi,E_chi] = histcounts(chi_map,hbins_chi);
            dE_chi = diff(E_chi(1:2));
            W_chi = V_chi/sum(V_chi);
            chi_mean = sum(E_chi(1:end-1).*W_chi);
            chi_std = sqrt(sum(W_chi.*(E_chi(1:end-1) + dE_chi/2 - chi_mean).^2));
            clim(ax3,[max(0,chi_mean-1*chi_std),chi_mean+1*chi_std])

            %%% Colorbar
            colormap(ax3_temp,'turbo')
            cb = colorbar(ax3_temp);
            cb.Title.String = '\chi^{2}';
            cb.FontSize = 15 ; 
            clim(ax3_temp,[max(0,chi_mean-1*chi_std),chi_mean+1*chi_std])
        else
            %%% Mask overlayed on intensity
            imagesc(ax2_2,Mask,'AlphaData', 0.4*Mask);
            title(ax2,"Peak Intensity",'FontSize',15)
            title(ax3,"No Lifetime Data Found",'FontSize',15)
            ax2_2.UserData = linkprop([ax2,ax2_2],{'x','y','Position','xlim'});

        end
        set(ax2,'XColor', 'none','YColor','none')
        set(ax2_2,'XColor', 'none','YColor','none')
        set(ax3,'XColor', 'none','YColor','none')
        set(ax3_temp ,'XColor', 'none','YColor','none')
        ax3_temp.Visible = 'off';
        ax2.Visible = 'off';
        ax2_2.Visible = 'off';
        ax3_temp.Visible = 'off';
        ax2_temp.Visible = 'off';

        Printing_Mask = double(Mask);
        [Img_Big_Seg, N_objects] = bwlabel(Mask,8);
        Obj_ptr = 1;
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
        xlim(ax,[xmin,xmax]);

        uicontrol('Parent',fig,'Style','text','String',"Label",'Position',[30 680 100 20]);
        listbox_ptr_label = uicontrol('Parent',fig,'Style','listbox','String',["Delete","1"],'Position',[30 380 100 300]);
        uicontrol('Parent',fig,'Style','text','String',"Type",'Position',[30 300 100 20]);
        listbox_ptr_type = uicontrol('Parent',fig,'Style','listbox','String',["Soma","Axon","Background","Cell 1","Cell 2"],'Position',[30 200 100 100]);
        btn_back = uicontrol('Parent',fig,'Style', 'pushbutton','String',"Back"  ,'Position', [200 10 100 35],'Callback',@FnBackBtn,'Enable','off');
        uicontrol('Parent',fig,'Style', 'pushbutton','String',"Next"  ,'Position', [300 10 100 35],'Callback',{@FnNextBtn,0});
        uicontrol('Parent',fig,'Style', 'pushbutton','String',"Delete"  ,'Position', [400 10 100 35],'Callback',{@FnNextBtn,-1});
        % uicontrol('Parent',fig,'Style', 'pushbutton','String',"Back"  ,'Position', [400 10 100 35],'Callback',{@FnBackBtn,-1});
        uicontrol('Parent',fig,'Style', 'pushbutton','String',"1"  ,'Position', [525 10 20 35],'Callback',{@FnNextBtn,1});
        uicontrol('Parent',fig,'Style', 'pushbutton','String',"Add Label"  ,'Position', [30 330 100 35],'Callback',@AddLabel);
        listbox_ptr_label.Value = 2;
        uiwait(fig)

        function AddLabel(~,~)
            lst_str = string(listbox_ptr_label.String);
            N = length(lst_str);
            if (N<20)
                listbox_ptr_label.Value = 1:N;
                listbox_ptr_label.String = ["Delete",string(1:N)];
                str_list = [str_list;strcat("Label ",num2str(N))];
                N_labels = length(str_list);
                Mycolormap_temp  = [Mycolormap;colorarr(1:N,:)];
                colormap(ax,Mycolormap_temp);
                colorbar(ax,'Ticks',1:N_labels,'Color','Black','FontSize',17,'TickLabels',str_list);
                clim(ax,[0,N_labels]);
                uicontrol('Parent',fig,'Style', 'pushbutton','String',num2str(N)  ,'Position', [500+N*25 10 25 35],'Callback',{@FnNextBtn,N});
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

        function FnNextBtn(~,~,var)
            if (Obj_ptr == 1)
                btn_back.Enable  = 'on';
            end
            Current_Type = listbox_ptr_type.Value;
            Current_Mask = (Img_Big_Seg == Obj_ptr);
            if (var == 0)   %% Next button was pressed
                Current_Label = listbox_ptr_label.Value-1;
            elseif(var == -1)
                Current_Label = 0;
            else
                Current_Label = var;
            end

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
            N_Obj = listbox_ptr_label.Value(end);
            delete(fig)
        end

    end
end


function [segment_Obj,Skip] = ThresholdingSlider_Multi(param, path_cell, Auto_threshold, folder_name_all, selected_filenames,filename_exist_arr)
% ThresholdingSlider_Multi
% Load peak-frame images for multiple measurement folders, then either:
%   - auto-threshold each independently (if Auto_threshold ~= 0), or
%   - open a per-image multi-control UI (sliderImageDemo_Multi) to set masks.
%
% Inputs:
%   param               : struct with fields like idx_peak_thresh (optional)
%   path_cell           : 1xN or Nx1 cell; each is a folder containing *.tif and RecSettings.txt
%   Auto_threshold      : 0 -> manual UI; otherwise numeric fraction of peak for auto mask
%   folder_name_all     : (unused for logic) pass-through for title context; we derive names from paths
%   selected_filenames  : 1xN cell of measurement folder names (for titles)
%
% Outputs:
%   segment_Obj         : 1xN cell of logical masks
%   Selected_ROI        : 1xN cell (empty; ROI not used here)

    % Normalize inputs
    if isrow(path_cell), path_cell = path_cell(:); end
    N_folder = numel(path_cell);

    if (param.Stain)
        path_main = fileparts(fileparts(path_cell{1}));
        if (param.current_stain_exists_flag)
            I_stain = [];
            try
                N_files = length(param.current_stain_filenames);
                for idx_stain = 1:N_files
                    Image_handle = Tiff(fullfile(path_main,param.current_stain_filenames{idx_stain}),'r');  % read given frame
                    I_stain(:,:,idx_stain) = read(Image_handle);
                end
            end
            if (~isempty(I_stain))
                [nx, ny, nz] = size(I_stain);

                % Reshape into 2D: each column is one image
                I_reshaped = reshape(I_stain, nx*ny, nz);
                
                % Convert to double (important for correlation)
                I_reshaped = double(I_reshaped);
                
                % Compute correlation matrix
                C = corrcoef(I_reshaped);
                
                % Apply threshold
                threshold = 0.9;
                is_good = C > threshold;
                mask = ~eye(nz); % true everywhere except diagonal
                % mean_corr = sum(C .* mask, 2) ./ sum(mask, 2);
                % bad_idx = find(mean_corr < threshold);
                % I_stain(:,:,bad_idx) = [];
                I_stain = mean(I_stain,3);
            end
        end
    end
    % Preallocate
    I = cell(N_folder,1);
    segment_Obj = cell(N_folder,1);
    Selected_ROI = cell(1, N_folder);  % keep shape consistent with other calls
    [Selected_ROI{:}] = deal([]);

    % ----- Load peak images exactly like ThresholdingSlider -----
    for idx_folder = 1:N_folder
        p = char(path_cell{idx_folder});
        files = dir([p, '/*.tif']);
        if isempty(files)
            I{idx_folder} = [];
            continue;
        end

        parametersFile = fullfile(p,'RecSettings.txt');
        [~,delta_t,tmin,tmax,~,BINNING,~,~] = GetExperimentalParameters(parametersFile); %#ok<ASGLU>
        t_sig = (tmin:delta_t:tmax)'; %#ok<NASGU>
        n = length(t_sig);

        % Sorted image list (same approach as your single-file loader)
        image_filenames = sort(split(strip(ls([p, '/*.tif']))));
        image_filenames = image_filenames(~cellfun('isempty',image_filenames)); % guard empties

        % Read peak frame (auto or fixed index)
        if isempty(param.idx_peak_thresh)
            I_temp = zeros(1216/BINNING, 1936/BINNING, n);
            for i = 1:n
                fn = erase(string(image_filenames{i}),"'");
                th = Tiff(fn,'r');
                I_temp(:,:,i) = double(read(th));
            end
            I_sum = squeeze(sum(I_temp,[1,2]));
            [~,idx_peak] = max(I_sum);
            I{idx_folder} = I_temp(:,:,idx_peak);
        else
            idxp = param.idx_peak_thresh;
            fn = erase(string(image_filenames{idxp}),"'");
            th = Tiff(fn,'r');
            I{idx_folder} = double(read(th));
        end
    end

    % ----- Path 1: Auto-threshold each independently -----
    if Auto_threshold ~= 0
        Tval = Auto_threshold;
        for idx_folder = 1:N_folder
            Img = imsharpen(I{idx_folder}, 'Radius', 30, 'Amount', 2);
            peak = max(Img(:));
            Mask_temp = Img > (Tval * peak);
            [Mask_temp, ~] = bwlabel(Mask_temp, 8);
            % Intensity floor for robustness (match your single version)
            Mask_temp(I{idx_folder} < 400) = 0;
            segment_Obj{idx_folder} = bwareaopen(Mask_temp, 50, 8) > 0;
        end
        Skip = 0;
        return;
    end

    % ----- Path 2: Manual UI with per-image controls (registration, CLim, etc.) -----
    % Use the dedicated multi-image UI you already have.
    % Titles: use derived day_names and your selected_filenames.
    [segment_Obj,Skip] = sliderImageDemo_Multi(I, folder_name_all, selected_filenames, filename_exist_arr,I_stain,param);
    if (param.Mask_together && Skip == 0)
        Mask_all = [];
        for idx_folder = 1:N_folder
            N_pixel = sum(sum(segment_Obj{idx_folder}));
            if (N_pixel > 0)
                if (isempty(Mask_all))
                    Mask_all = segment_Obj{idx_folder};
                else
                    Mask_all = Mask_all.*segment_Obj{idx_folder};
                end
            end
        end
        for idx_folder = 1:N_folder
            N_pixel = sum(sum(segment_Obj{idx_folder}));
            if (N_pixel > 0)
                segment_Obj{idx_folder} = Mask_all;
            end
        end
    end

end


function [Mask,Skip] = sliderImageDemo_Multi(Peak_Imgs, folder_name, selected_filenames, filename_exist_arr,I_stain,param)
% Multi-image thresholding + per-image registration UI (adapted from legacy code).
% Returns: 
%   Mask         : 1xN cell of logical masks (per image, with applied translation)
%   Selected_ROI : 1xN cell (empty; ROI not used here)
    Skip = 0;
    % -------- Normalize inputs to legacy shapes --------
        N_days = sum(filename_exist_arr);
        days_arr = find(filename_exist_arr == 1);
        if iscell(Peak_Imgs)
            for idx_days = 1:N_days
                idx_current = days_arr(idx_days);
                if (idx_days == 1)
                    sz1 = size(Peak_Imgs{idx_current},1); sz2 = size(Peak_Imgs{idx_current},2);
                    I = zeros(sz1, sz2, N_days);
                end                
                I(:,:,idx_days) = double(Peak_Imgs{idx_current});
                folder_name_arr(idx_days) = string(folder_name{idx_current});
                selected_filenames_arr(idx_days) = string(selected_filenames{idx_current});
            end
        else
            I(:,:) = double(Peak_Imgs); 
            folder_name_arr = string(folder_name);
            selected_filenames = string(selected_filenames);
        end
        
    % -------- Compute truncation mask like the old code --------
    refIdx   = min(2, N_days);
    Img_temp = I(:,:,refIdx);
    Img_trunc = (Img_temp < 200);
    Img_trunc = bwmorph(Img_trunc,'bridge',inf);
    Img_trunc = bwareaopen(Img_trunc,50,4);
    Img_trunc = 1-Img_trunc;
    Img_trunc = bwareaopen(Img_trunc,50,4);

    % -------- UI state (names & defaults mirror legacy) --------
    fig = uifigure;                 
    flag_show_mixture = 0;
    if (flag_show_mixture)
        gcf7 = figure(7); %#ok<LFIG>
    end
    ax   = cell(N_days,1);    ax2  = cell(N_days,1);
    sliders        = cell(N_days,1);
    sliders_filter = cell(N_days,1);
    sliders_clim_min = cell(N_days,1);
    sliders_clim_max = cell(N_days,1);
    label_clim   = cell(N_days,1);
    label        = cell(N_days,1);
    label_filter = cell(N_days,1);
    T   = cell(N_days,1);
    N_f = cell(N_days,1);
    imgHandle0   = cell(N_days,1);
    imgHandleStain = cell(N_days,1);
    imgHandleR   = cell(N_days,1);
    ax_stain = [];
    imgHandle_stain = [];
    sliders_x    = cell(N_days,1);
    sliders_y    = cell(N_days,1);
    checkbox     = cell(N_days,1);
    rb           = cell(N_days,1);
    X_pos        = zeros(N_days,1);
    Y_pos        = zeros(N_days,1);
    Used_mask        = cell(N_days,1);
    Used_mask_temp   = cell(N_days,1);

    % Display & control defaults (as in legacy)
    Thresh   = 0.15;
    N_filter = 100;
    h = 10;  w = 50;
    h2 = -10; w2 = 100;
    x0=10; y0=200; width=960; height=600;

    % Normalize images to [0,1] for display
    for i = 1:N_days
        Ii = I(:,:,i);
        I(:,:,i) = Ii;
    end

    % -------- Layout (kept very close to legacy) --------
    fig.Position  = [50, 50, 1700, 800];
    frame1_new  = uipanel(fig,'Position',[50, 290, 1500, 120]);
    frame2_new  = uipanel(fig,'Position',[50, 20, 1500, 250]);
    fig.AutoResizeChildren  = 'off';
    frame2_new.AutoResizeChildren = 'off';
    if (flag_show_mixture)
        set(gcf7,'Color','k'); set(gcf7,'position',[x0,y0,width,height]);
    end

    uicontrol('Parent',frame1_new,'Style','text','Position',[w, h+60, 100, 20],'String','Threshold');
    uicontrol('Parent',frame1_new,'Style','text','Position',[w, h+30, 100, 20],'String','Filter');
    uicontrol('Parent',frame1_new,'Style','text','Position',[w, h,    100, 20],'String','CLim');

    bg = uibuttongroup(frame2_new,'Position',[20 80 123 130],'Title','Reference');
    dx1_slider = 260; dx2_slider = 210; dx_img = 0.25;

    % One hidden axes per image in gcf7 for colored masks + colorbar linkage
    uicontrol('Parent',fig,'Style','text','Position',[400, 775, 1000, 25],'String',selected_filenames{1},'FontSize',15);
    for i = 1:N_days
        % Controls per image
        sliders{i}        = uicontrol('Parent',frame1_new,'Style','slider','Min',0,'Max',1,'Value',Thresh, ...
                                      'Position',[w+100+(i-1)*dx1_slider, h+60, 110, 20], 'Callback',@updateImage);
        sliders_filter{i} = uicontrol('Parent',frame1_new,'Style','slider','Min',0,'Max',120,'Value',N_filter, ...
                                      'Position',[w+100+(i-1)*dx1_slider, h+30, 110, 20], 'Callback',@updateImage);
        sliders_clim_min{i} = uicontrol('Parent',frame1_new,'Style','slider','Min',200,'Max',500,'Value',220, ...
                                      'Position',[w+100+(i-1)*dx1_slider, h, 55, 20], 'Callback',@updateImage);
        sliders_clim_max{i} = uicontrol('Parent',frame1_new,'Style','slider','Min',1000,'Max',3000,'Value',1000, ...
                                      'Position',[w+155+(i-1)*dx1_slider, h, 55, 20], 'Callback',@updateImage);
        label_clim{i}     = uicontrol('Parent',frame1_new,'Style','text','Position',[w+220+(i-1)*dx1_slider, h, 100, 20], ...
                                      'String','[220,1000]');
        label{i}          = uicontrol('Parent',frame1_new,'Style','text','Position',[w+220+(i-1)*dx1_slider, h+60, 60, 20], ...
                                      'String',num2str(Thresh));
        label_filter{i}   = uicontrol('Parent',frame1_new,'Style','text','Position',[w+220+(i-1)*dx1_slider, h+30, 60, 20], ...
                                      'String',num2str(N_filter));
        T{i}   = Thresh;
        N_f{i} = N_filter;

        % Overlay axes on fig (intensity + mask overlay)
        ax2{i} = axes(fig);
        ax2{i}.Units = 'pixels';
        ax2{i}.Position = [30+(i-1)*310, 450, 300, 300];
        Im_plot{i} = (min(max(I(:,:,i),220),1000)-220)./(1000-220);
        Img_overlay = labeloverlay(Im_plot{i}, i*Thresh_Image(Thresh, I(:,:,i), N_filter, Img_trunc), 'Colormap', jet(N_days),'Transparency',0);
        imgHandle0{i} = imshow(Img_overlay,'Parent',ax2{i});
        if (param.Stain)
            if (param.current_stain_exists_flag)
                hold(ax2{i}, 'on')
                St = double(I_stain);
                St = St ./ max(St(:));
                stainRGB = cat(3, ones(size(St)), zeros(size(St)), zeros(size(St)));
                imgHandleStain{i} = imshow(stainRGB, 'Parent', ax2{i});
                imgHandleStain{i}.AlphaData = 0.5 * St;
            end
        end
        xlim(ax2{i},[90,770]); axis(ax2{i},'square');
        set(ax2{i},'XColor','none','YColor','none');
        colormap(ax2{i}, gray(256));
        title(ax2{i}, sprintf('File %d: %s', i, folder_name_arr(i)), 'Interpreter','none');

        % Colored mask axes on gcf7 (hidden but used for mask layer/colorbar)
        if (flag_show_mixture)
            if i==1
                ax{i} = axes(gcf7);
                ax{i}.Position = [0.1300 0.1100 0.7750 0.8150];
            else
                ax{i} = copyobj(ax{1}, gcf7);
                ax{i}.Position = [0.1300 0.1100 0.7750 0.8150];
                ax{i}.UserData = linkprop([ax{1}, ax{i}], {'x','y','Position'});
            end
            set(ax{i},'XColor','none','YColor','none','Visible','off');
        end
        Mk0 = Thresh_Image(Thresh, I(:,:,i), N_filter, Img_trunc);
        if (flag_show_mixture)
            imgHandleR{i} = imagesc(ax{i}, i*Mk0, 'AlphaData', Mk0);
            set(ax{i},'XColor','none','YColor','none','Visible','off');
            clim(ax{i},[1,N_days]); colormap(ax{i}, jet(N_days));
        end
    
        % colorbar(ax{i}, 'Ticks',1:N_days,'Color','White','FontSize',17, 'TickLabels',string(1:N_days));
        % Registration controls
        sliders_x{i} = uislider('Parent',frame2_new,'Limits',[-400 400],'Value',0, ...
                                 'Position',[w2+110+(i-1)*dx2_slider, h2+60, 150, 20], ...
                                 'ValueChangingFcn',@updateImageReg,'ValueChangedFcn',@updateImageReg);
        sliders_y{i} = uislider('Parent',frame2_new,'Orientation','Vertical','Limits',[-400 400],'Value',0, ...
                                 'Position',[w2+100+(i-1)*dx2_slider, h2+80, 20, 150], ...
                                 'ValueChangingFcn',@updateImageReg,'ValueChangedFcn',@updateImageReg);
        checkbox{i}  = uicheckbox(frame2_new,'Text','Show Mask', ...
                                  'Position',[w2+180+(i-1)*dx2_slider, h2+150, 100, 20], ...
                                  'ValueChangedFcn',@toggleMask,'Value',1);
        rb{i}        = uiradiobutton(bg,'Position',[10, 100-20*i, 100, 20],'Text',sprintf('Day %d',i));

        % Initialize masks
        Used_mask{i}      = Mk0;
        Used_mask_temp{i} = Mk0;
    end

    if (param.Stain)
        if (param.current_stain_exists_flag)
            I_stain_plot = double(I_stain);
            I_stain_plot = I_stain_plot - min(I_stain_plot(:));
            I_stain_plot = I_stain_plot ./ max(I_stain_plot(:));
            
            ax_stain = axes(fig);
            ax_stain.Units = 'pixels';
            ax_stain.Position = [30 + N_days*310, 450, 300, 300];
            
            imgHandle_stain = imshow(I_stain_plot, 'Parent', ax_stain);
            axis(ax_stain,'square');
            set(ax_stain,'XColor','none','YColor','none');
            title(ax_stain, 'HUC stain');
            xlim(ax_stain,[90,770]);
        end
    end
    sliders_x{1}.Enable = 0; sliders_y{1}.Enable = 0; rb{1}.Value = true;

    % Action buttons
    uicontrol('Parent',frame1_new,'Style','pushbutton','String','Update Reg', ...
              'Position',[w+350+(N_days-1)*dx1_slider, h+55, 100, 35],'Callback',@SetReg);
    uicontrol('Parent',frame1_new,'Style','pushbutton','String','Skip', ...
          'Position',[w+350+(N_days-1)*dx1_slider, h+15, 100, 35],'Callback',@SkipThis);
    uicontrol('Parent',frame2_new,'Style','pushbutton','String','Start Analyzing', ...
              'Position',[w+350+(N_days-1)*dx2_slider, 70, 100, 35],'Callback',@SetT);
    uicontrol('Parent',frame2_new,'Style','pushbutton','String','Auto Register', ...
              'Position',[w+350+(N_days-1)*dx2_slider, 120, 100, 35],'Callback',@SetRegAuto);
    uicontrol('Parent',frame2_new,'Style','pushbutton','String','Change Reference', ...
              'Position',[30, 15, 100, 35],'Callback',@SetRef);

    uiwait(fig);  % block until Start Analyzing

    % -------- Nested fns (names & logic closely track legacy) --------
    function mask_temp = Thresh_Image(Tin, Peak_Img, Nfilter, ImgTr)
        Peak_ImgS  = imsharpen(Peak_Img,'Radius',30,'Amount',2);
        Peak_ImgS  = adapthisteq(ImgTr .* Peak_ImgS ./ max(1, max(Peak_ImgS(:))));
        mask_temp  = imbinarize(Peak_ImgS, Tin);
        mask_temp  = bwareaopen(mask_temp, Nfilter, 8);
        mask_temp  = bwmorph(mask_temp,'bridge',inf);
        mask_temp  = bwmorph(mask_temp,'fill',1);
        mask_temp(Peak_Img < 100) = 0;  % robustness as in your newer code
    end

    function updateImage(~,~)
        for i = 1:N_days
            Thresh   = get(sliders{i}, 'Value');
            N_filter = round(get(sliders_filter{i}, 'Value'));
            minValue = get(sliders_clim_min{i}, 'Value');
            maxValue = get(sliders_clim_max{i}, 'Value');

            if (T{i} ~= Thresh) || (N_f{i} ~= N_filter) || ...
               (abs(minValue - getVal(label_clim{i},1))>0 || abs(maxValue - getVal(label_clim{i},2))>0)
                T{i}   = Thresh;   N_f{i} = N_filter;
                Mk     = Thresh_Image(Thresh, I(:,:,i), N_filter, Img_trunc);
                % Update overlay with CLim scaling
                Scaled_I = (min(max(I(:,:,i),minValue),maxValue)-minValue)./(maxValue-minValue);
                Img_overlay = labeloverlay(Scaled_I, i*Mk, 'Colormap', jet(N_days),'Transparency',0);
                set(imgHandle0{i}, 'CData', Img_overlay);
                set(label_clim{i}, 'String', sprintf('[%d,%d]', round(minValue), round(maxValue)));
                set(label{i}, 'String', sprintf('%.3f', Thresh));
                set(label_filter{i}, 'String', num2str(N_filter));
            end
        end
    end

    function v = getVal(lbl, idx)
        % parse "[a,b]" string safely
        s = string(get(lbl,'String'));
        s = erase(s,{'[',']'});
        parts = split(s,',');
        if numel(parts) < 2, v = 0; return; end
        v = str2double(parts(idx));
        if isnan(v), v = 0; end
    end

    function updateImageReg(~,~)
        for i = 1:N_days
            x_axis = round(get(sliders_x{i}, 'Value'));
            y_axis = round(get(sliders_y{i}, 'Value'));
            if (X_pos(i) ~= x_axis) || (Y_pos(i) ~= y_axis)
                X_pos(i) = x_axis; Y_pos(i) = y_axis;
                Used_mask_temp{i} = imtranslate(Used_mask{i}, [x_axis, y_axis]);
                set(imgHandleR{i}, 'CData', i*Used_mask_temp{i});
                set(imgHandleR{i}, 'AlphaData', Used_mask_temp{i});
                if (flag_show_mixture)
                    colormap(ax{i}, jet(N_days));
                end
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

    function SetReg(~,~)
        for i = 1:N_days
            Thresh   = get(sliders{i}, 'Value');
            N_filter = round(get(sliders_filter{i}, 'Value'));
            Used_mask{i} = Thresh_Image(Thresh, I(:,:,i), N_filter, Img_trunc);
            set(imgHandleR{i}, 'CData', i*Used_mask{i});
            set(imgHandleR{i}, 'AlphaData', Used_mask{i});
        end
    end

    function SetRegAuto(~,~)
        % Find current reference
        Ref = 1;
        for i = 1:N_days
            if rb{i}.Value, Ref = i; break; end
        end
        for i = 1:N_days
            if i == Ref, continue; end
            tform  = imregcorr(double(Used_mask_temp{i}), double(Used_mask_temp{Ref}), "translation");
            x_axis = tform.Translation(1);
            y_axis = tform.Translation(2);
            Used_mask_temp{i} = logical(imtranslate(Used_mask{i}, [x_axis, y_axis]));
            set(imgHandleR{i}, 'CData', i*Used_mask_temp{i}, 'AlphaData', Used_mask_temp{i});
            sliders_x{i}.Value = x_axis; sliders_y{i}.Value = y_axis;
        end
    end

    function SetRef(~,~)
        for i = 1:N_days
            if rb{i}.Value
                Used_mask_temp{i}   = Used_mask{i};
                sliders_x{i}.Enable = 0; sliders_y{i}.Enable = 0;
                sliders_x{i}.Value  = 0; sliders_y{i}.Value  = 0;
                if checkbox{i}.Value
                    set(imgHandleR{i}, 'AlphaData', Used_mask_temp{i}, 'CData', i*Used_mask_temp{i});
                end
            else
                sliders_x{i}.Enable = 1; sliders_y{i}.Enable = 1;
            end
        end
    end

    function SetT(~,~)
        % Finalize and return per-image masks with current translation
        for i = 1:N_days
            Thresh   = get(sliders{i}, 'Value');
            N_filter = round(get(sliders_filter{i}, 'Value'));
            Mk = Thresh_Image(Thresh, I(:,:,i), N_filter, Img_trunc);
            if rb{i}.Value
                x_axis = 0; y_axis = 0;
            else
                x_axis = round(get(sliders_x{i}, 'Value'));
                y_axis = round(get(sliders_y{i}, 'Value'));
            end
            Mk = imtranslate(Mk, [x_axis, y_axis]);
            Mask{i} = logical(Mk);
        end
        delete(fig); delete(frame2_new); 
        if (flag_show_mixture)
            if ishghandle(gcf7), close(gcf7); end
        end
        % uiresume; % safety
    end

    function SkipThis(~,~)
        % Finalize and return per-image masks with current translation
        Mask = [];
        Skip = 1;
        delete(fig); delete(frame2_new);
        if (flag_show_mixture)
            if ishghandle(gcf7), close(gcf7); end
        end
        % uiresume; % safety
    end
end

%% Fitting
function [PriorCell,Returned_data,status] = SalemPixelFitting22(arg1, arg2, param)

    path = char(arg1);
    status = 1;
    Returned_data = [];
    PriorCell = [];
    %%% Setting file name
    N_files = size(path,1);
    I = cell(N_files,1);
    err_map = cell(N_files,1);
    if (param.IRF_1D == 0)
        IRF_file = strcat('IRFData/IRF',num2str(param.IRF_dt),'t_',num2str(200),'G_',num2str(2),'BIN_',param.laser,'.mat');
    else
        IRF_file = strcat('IRFData/IRF',num2str(param.IRF_dt),'t_',num2str(200),'G_',num2str(2),'BIN_',param.laser,'_1D.mat');
    end
    if isfile(IRF_file)
        load(IRF_file,'IRF','time_axis_IRF','t_IRF');   
        if (exist('time_axis_IRF','var'))
            t_IRF = time_axis_IRF';
        else
            t_IRF = t_IRF';
        end
    else
        disp('IRF Does not exist!')
        disp(strcat('Expected file: ',IRF_file, ', in path: ',pwd))
        status = 3;
        return
    end

    filename_saving_results = Get_savefile_name(param,path(1,:),"Results");

    if (param.Stain == 1)
        if (param.current_stain_exists_flag)
            try
                N_files = length(param.current_stain_filenames);
                for idx_stain = 1:N_files
                    Image_handle = Tiff(fullfile(fileparts(fileparts(path(1,:))),param.current_stain_filenames{idx_stain}),'r');  % read given frame
                    I_stain(:,:,idx_stain) = read(Image_handle);
                end
                param.stain_success = 1;
                I_stain = mean(I_stain,3);
            catch
                param.stain_success = 0;
            end
        end
    end

    for idx_path = 1:size(path,1)
        if (path(idx_path,end) == '/')
            path(idx_path,:) = path(idx_path,1:end-1);
        end
        parametersFile{idx_path} = strcat(path(idx_path,:),'/RecSettings.txt'); % Contains experimental parameters
        param_cpd{idx_path} = param;
        [~,delta_t,tmin,tmax,~,BINNING,MCP,SubDrk] = GetExperimentalParameters(parametersFile{idx_path});    % get experimental parameters
        Noise_file = strcat('NoiseData/Noise_histogram_',num2str(MCP),'_v',num2str(param.NoiseVersion),'.mat');
        if exist(Noise_file,'file')
            load(Noise_file,param.NoiseFileName)
        else
            fprintf("Noise file not found: %s",Noise_file);
            status = 27;
            return
        end
        param_cpd{idx_path}.SubDrk = SubDrk;
        param_cpd{idx_path}.MCP = MCP;
        param_cpd{idx_path}.MCPNoise = eval(param.NoiseFileName);
        HorizontalPixels = 1:1936/BINNING;     % based on binning factor [pixels]
        VerticalPixels = 1:1216/BINNING;         % based on binning factor [pixels]
        t_sig{idx_path} = (tmin:delta_t:tmax)';
        n = length(t_sig{idx_path});
        image_filenames = sort(split(strip(ls([path(idx_path,:), '/*.tif']))));   % image files (sorted)
        I{idx_path} = zeros(range(VerticalPixels)+1,range(HorizontalPixels)+1,n);       % matrix which array C will be converted to
        I_temp = zeros(range(VerticalPixels)+1,range(HorizontalPixels)+1,n);
        try
            for i = 1:n     %% extract pixel intensities and time delays for each frame
                image_filename_temp = erase(string(image_filenames{i}) , "'" );
                Image_handle = Tiff(image_filename_temp,'r');  % read given frame
                I_temp(:,:,i) = read(Image_handle); % extract intensity data for given frame
            end
            I{idx_path} = I_temp;
        catch
            status = 0;
            disp('Cannot load Tif file!')
            return;     %%% End the program if something is wrong in the file
        end
        I_sum = sum(I_temp,[1,2]);
        [~,idx_peak] = max(I_sum);
        Volt_gain_relation  = [0.0054	0.00755	0.0097	0.01185	0.014	0.0206	0.0272	0.0338	0.0404	0.047	0.0656	0.0842	0.1028	0.1214	0.14	0.191	0.242	0.293	0.344	0.395	0.5214	0.6478	0.7742	0.9006	1.027	1.3336	1.6402	1.9468	2.2534	2.56	3.254	3.948	4.642	5.336	6.03	7.578	9.126	10.674	12.222	13.77	16.856	19.942	23.028	26.114	29.2	34.9	40.6	46.3	52	57.7	68.275	78.85	89.425	100];
        Volt_axis = 260:10:790;
        MCP_idx = (Volt_axis==MCP);
        gain_ratio = Volt_gain_relation(MCP_idx);
        param_cpd{idx_path}.idx_peak = idx_peak;
        param_cpd{idx_path}.gain_ratio = gain_ratio;
        [~,idx] = max(sum(I{idx_path},[1 2]));
        
        if (N_files > 1)
            param_cpd{idx_path}.Mask = arg2{idx_path};
            Peak_Img{idx_path} = I{idx_path}(:,:,idx);
            I_avg{idx_path} = sum(I{idx_path},3)/size(I{idx_path},3);
        else
            Peak_Img = I{idx_path}(:,:,idx);
            I_avg = sum(I{idx_path},3)/size(I{idx_path},3);
        end
        err_map{idx_path} = zeros(size(I{idx_path},1),size(I{idx_path},2)); %% Err Map is the main vector that holds the information about which pixels to fit or not
        err_map{idx_path}(param.Mask == 0) = 1;
        [status] = CheckTimeAxis(t_sig{idx_path},t_IRF);
        if (status == 2)
            disp("Error in time axis")
            return
        end
        if (sum(param_cpd{idx_path}.Mask(:)) == 0)
            disp('Not enough pixels. Check thresholding!')
            status = 4;
            return
        end
        %%%% Jump
        Obj_Mask_temp = [];
        if (param.Analysis_type == 3 && param.Stain == 1)
            if (param.stain_success)
                %%%%%
                [Stain_mask,Size_mask,cellsize_arr,cellsize_arr,Mask] = getHealthTypeAndSize(I_stain,param_cpd{idx_path}.Mask);
                Obj_Mask_temp = Stain_mask;
            end
        end

        if(isempty(Obj_Mask_temp))
            Mask_common_filtered = bwareaopen(param.Mask ,50,8);
            [Obj_Mask,N_Object] = bwlabel(Mask_common_filtered,8);
        else
            Obj_Mask = Obj_Mask_temp;
        end


        if(param.AmpEst_flag == 1)
            AmpMap_path = Get_savefile_name(param,path,"AmpMap");
            load(AmpMap_path,"Amp_Map")
            param.AmpMap = Amp_Map;
        end
    end
    [status,Returned_data] = TimeDomainFit(I,IRF,t_sig,t_IRF,err_map,Obj_Mask,param_cpd);
    if (status >= 5)
        if (status == 17)
            disp('Debugging complete')
        else
            disp('Errors during the fitting')
        end
        return
    end
    Returned_data.time = t_sig;
    Returned_data.I_avg = I_avg;
    Returned_data.Peak_Img = Peak_Img;
    Returned_data.HorizontalPixels = HorizontalPixels;
    Returned_data.VerticalPixels = VerticalPixels;
    Returned_data.Mask_used = param.Mask;
    Returned_data.Photon_Budget = sum(I_temp-param.DCShift,3)/gain_ratio;
    Returned_data.MCP = MCP;
    Returned_data.Total_counts = sum(I_temp-param.DCShift,3);
    Returned_data.Gain = gain_ratio;

    if (param.Stain)
        if (param.current_stain_exists_flag)
            if (param.stain_success)
                Returned_data.StainImg = I_stain;
                if (param.Analysis_type == 3)
                    Returned_data.Size_mask = Size_mask;
                end
            end
        end
    end
    save(filename_saving_results,'Returned_data','-mat')
    if isfield(Returned_data,'tau_map')
        PriorCell = Returned_data.tau_map;
    else
        PriorCell = [];
    end
end
%% Functions for processing the data
function [status] = CheckTimeAxis(t_sig,t_IRF)
%%% This will only work correctly if I range is includede in IRF range
    status = 1;    
    dt_IRF = t_IRF(2)-t_IRF(1);
    dt = t_sig(2)-t_sig(1);
    
    if (dt < dt_IRF && t_IRF(1) > t_sig(1) && t_IRF(end) < t_sig(end))
        fprintf("Bad time axis, Signal needs to be clipped!")
        status = 2;
    elseif (t_IRF(1) > t_sig(1) && t_IRF(end) < t_sig(end))
        fprintf("Bad time axis, signal will be clipped!")
        status = 2;
    end    
end

function [status,Returned_data] = TimeDomainFit(I,IRF, t_sig, t_IRF, err_map ,Label_Mask, param_folder)
    status = 1;    %%%% No error until otherwise stated
    err_str = [];
    Noise_calib = [];
    Noise_Var_vec = [];
    Prior_Current = [];
    Returned_data = struct();
    Returned_data.flag_empty   = 0;
    param = param_folder{1};
    N_path = size(param_folder,1);
    [Nx,Ny,Nt] = size(I{1});
    if (param.pixelwise == 1)
        amp_map = zeros(Nx,Ny,param.order);       
        chi_map = zeros(Nx,Ny);
        Ar_map = zeros(Nx,Ny);
        DC_map = zeros(Nx,Ny);
        tshift_map = zeros(Nx,Ny);
        tau_map = zeros(Nx,Ny,param.order);
    else
        chi_map = [];
        tau_map = [];
        amp_map = [];
    end

    if(~isempty(param.ROI)), ROI_map = zeros(Nx,Ny); ROI_map(param.ROI(2):param.ROI(2)+param.ROI(4),param.ROI(1):param.ROI(1)+param.ROI(3)) = 1; end

    if (param.IRF_1D == 0), IRF_current = squeeze(sum(IRF,[1 2]));
    else, IRF_current = IRF; end

    if (param.NeighborsNoise_flag == 1 || param.Method == 7)    
        %%% Salem: I just moved things but never tested this.
        Noise_calib = cell(N_path,1);
        for idx_folder = 1:N_path
            Noise_calib{idx_folder} = Get_Noise_Neighbors(param,I{idx_folder});
        end
    end
    if (param.pixelwise == 1)    %%% Fitting Pixles
        %%%%%%%%%%%%
        %%% We assume single file from here on since pixel-wise is a single
        %%% file operation
        if (~isempty(param.ROI)), [rows, cols] = find(ROI_map > 0 & err_map{1} == 0);
        else, [rows, cols] = find(err_map{1} == 0); end
        idx_vector = [rows, cols];
        N_pixels = length(rows);
        fprintf("%d pixels. ",N_pixels);
        chi_map_vec = cell(N_pixels,1);
        Ar_map_vec = cell(N_pixels,1);
        amp_map_vec = cell(N_pixels,1);
        tau_map_vec = cell(N_pixels,1);
        err_vec = cell(N_pixels,1);
        err_str_vector = cell(N_pixels,1);
        err_status_vec = cell(N_pixels,1);
        I_fit   = cell(N_pixels,1);
        warning('off','all');
        warning('off','MATLAB:remoteparfor:ParforWorkerAborted');
        if (N_pixels == 0)
            disp('No points to analyze! Exiting!');
            status = 7;
            return;
        end
        count_bad = 0;

                %%%%% Form noise matrix
        x_indices = (1:Nx).';
        y_indices = 1:Ny;
            
        %%% If spatial neighbor distance is needed
        [rows, cols] = find(param.Mask > 0);
        N_pixels = length(rows);
        W_local = zeros(N_pixels, Nt);      % worker-safe sliced variable
        fprintf("Estimating Noise...")
        parfor idx_pixel = 1:N_pixels
            i = rows(idx_pixel);j = cols(idx_pixel);
            d = ((x_indices-i).^2+(y_indices-j).^2);
            [~,idx] = sort(d(:),'ascend');
            nn_lin = idx(1:param.N_neighbors);
            [i_near,j_near] = ind2sub([Nx,Ny],nn_lin);
            y_neighbors = I{1}(i_near,j_near,:);
            y_neighbors_avg = squeeze(mean(y_neighbors,[1,2]));
            y_neighbors_smoothed = smooth(y_neighbors_avg, 7);
            NoiseVector_empirical  = param.MCPNoise(cast(abs(y_neighbors_smoothed),"int16"));
            W_local(idx_pixel, :) = (1 ./ NoiseVector_empirical).';
        end
        W{1} = zeros(Nx,Ny,Nt);
        for p = 1:N_pixels
            W{1}(rows(p), cols(p), :) = reshape(W_local(p, :), 1, 1, Nt);
        end


        if (param.DEBUGGING == 0)
            tic
            for count = 1:min(N_pixels,10)
                i = rows(count);j = cols(count);
                y_current = squeeze(I{1}(i,j,:));
                w_current = squeeze(W{1}(i,j,:));
                x_indices = rows-i;y_indices = cols-j;
                d = (x_indices.^2+y_indices.^2);
                [~,idx] = sort(d);
                if (length(rows) >= param.N_neighbors)
                    nearest_x = rows(idx(2:param.N_neighbors+1));nearest_y = cols(idx(2:param.N_neighbors+1));
                else
                        
                    % if (length(rows) > 1)
                    %     nearest_x = rows(idx(2:end));nearest_y = cols(idx);
                    % else
                    %     nearest_x = rows(idx);nearest_y = cols(idx);
                    % end
                end
                linear_indices = sub2ind([Nx,Ny], nearest_x, nearest_y);
                I_temp = reshape(I{1},[Nx*Ny,Nt]);
                W_temp = reshape(W{1},[Nx*Ny,Nt]);
                y_neighbors = I_temp(linear_indices, :);
                w_neighbors = W_temp(linear_indices, :);
                nearest_x_tshift = rows(idx(1:min(N_pixels,20)));
                nearest_y_tshift = cols(idx(1:min(N_pixels,20)));
                linear_indices_tshift = sub2ind([Nx, Ny], nearest_x_tshift, nearest_y_tshift);
                y_avg_tshift = mean(I_temp(linear_indices_tshift, :));
                if ~isempty(param.PriorData)
                    Prior_temp = cell(param.order,1);
                    Prior_temp{2} = param.PriorData(i,j);
                    Prior_Current{1} = Prior_temp;
                end
                if (param.AmpEst_flag == 1)
                    Amp_temp = param.AmpMap(i,j);
                    Prior_Current{2} = [{1-Amp_temp},{Amp_temp}];
                end
                if (param.NeighborsNoise_flag == 1 || param.Method == 7), NoiseVector = squeeze(Noise_calib(i,j,:));
                else, NoiseVector = []; end
                try
                    mydeconv8(y_current, w_current, t_sig{1}, IRF_current, t_IRF, param, NoiseVector, y_neighbors, w_neighbors, y_avg_tshift, Prior_Current);
                catch
                    count_bad = count_bad + 1;
                end
            end
            if (count_bad > min(N_pixels,10))
                disp('Too Many Errors!')
                status =  5;
                return
            end
            elapsedTime = toc;
            ETA = (elapsedTime / count) * (N_pixels - count) / 60;
            try
                pool = gcp;
                numWorkers = pool.NumWorkers;
                fprintf('ETA: %.2f (%d), ', round(ETA/numWorkers,2),numWorkers);
            catch
                disp(strcat('ETA: ',num2str(round(ETA,2)),' mins'))
            end

            tic
            dq = parallel.pool.DataQueue;
            afterEach(dq, @updateProgress);
            fprintf('Progress:  ');
            parfor count = 1:N_pixels
                Prior_Current = [];
                i = rows(count);j = cols(count);
                err_str = [];
                y_current = squeeze(I{1}(i,j,:));      %%%% Future update: I can be stretched to a  2D vector of good pixels only to save memory
                w_current = squeeze(W{1}(i,j,:));
                x_indices = rows-i;y_indices = cols-j;
                d = (x_indices.^2+y_indices.^2);
                [~,idx] = sort(d);
                nearest_x = rows(idx(2:param.N_neighbors+1));nearest_y = cols(idx(2:param.N_neighbors+1));
                linear_indices = sub2ind([Nx, Ny], nearest_x, nearest_y);
                I_temp = reshape(I{1},[Nx*Ny,Nt]);
                W_temp = reshape(W{1},[Nx*Ny,Nt]);
                y_neighbors = I_temp(linear_indices, :);
                w_neighbors = I_temp(linear_indices, :);
                nearest_x_tshift = rows(idx(1:min(N_pixels,20)));nearest_y_tshift = cols(idx(1:min(N_pixels,20)));
                linear_indices_tshift = sub2ind([Nx, Ny], nearest_x_tshift, nearest_y_tshift);
                y_avg_tshift = mean(I_temp(linear_indices_tshift, :));
                if (param.AmpEst_flag == 1)
                    Amp_temp = param.AmpMap(i,j);
                    Prior_Current{2} = [{1-Amp_temp},{Amp_temp}];
                end
                if ~isempty(param.PriorData)
                    Prior_temp = cell(param.order,1);
                    Prior_temp{2} = param.PriorData(i,j);
                    Prior_Current{1} = Prior_temp;
                end
                if (param.NeighborsNoise_flag == 1 || param.Method == 7)
                    NoiseVector = squeeze(Noise_calib(i,j,:));
                else
                    NoiseVector = [];
                end

                try
                    [y_fit,Noise_var,amp,tau,chi,err_status,chi_vec,Ar_out,DC_out,tshift_out,err_str,~] = mydeconv8(y_current,w_current, t_sig{1}, IRF_current, t_IRF, param, NoiseVector, y_neighbors, w_neighbors, y_avg_tshift, Prior_Current);
                    if (err_status == 0)     %%% If we caught an error
                        I_fit{count} = y_fit;
                        amp_map_vec{count} = amp;
                        % Noise_Var_vec{count} = Noise_var;
                        tau_map_vec{count} = tau;
                        chi_map_vec{count}  = chi;
                        chi_vec_arr{count}  = chi_vec;
                        Ar_map_vec{count}  = Ar_out;
                        DC_map_vec{count}  = DC_out;
                        tshift_map_vec{count}  = tshift_out;
                        err_vec{count} = 0;
                        err_str_vector{count} = [];
                    end
                    err_status_vec{count} = err_status;
                        % Debugging aid: save bad RCOND matrices
                catch ME  %%% If a new unknown error happened
                    err_status = 1;
                    err_str_vector{count}  = ME.message;
                    err_status_vec{count} = err_status;
                    disp(strcat('Error in fitting! At count: ',num2str(count)))
                    if contains(ME.message, 'singular')
                        warning('⚠️ Singular matrix at pixel %d: %s', count, ME.message);
                    end

                end
                if (err_status > 0)     %%% If we caught an error
                    I_fit{count} = [];
                    Noise_Var_vec{count} = NaN;
                    amp_map_vec{count} = NaN;
                    tau_map_vec{count} = NaN;
                    chi_map_vec{count}  = NaN;
                    Ar_map_vec{count}  = NaN;
                    DC_map_vec{count}  = NaN;
                    tshift_map_vec{count}  = NaN;
                    chi_vec_arr{count}  = [];
                    err_vec{count} = err_status;
                    err_str_vector{count}  = err_str;
                end
                err_status_vec{count} = err_status;
                send(dq, i); 
            end
        end

        if (param.DEBUGGING == 1)
            Debug_counter = 1;        %%% Insert a specific pixel number here for quick debugging
            for count = 1:N_pixels
                if (param.DEBUG_Count ~= 0)
                    count = param.DEBUG_Count;
                    if (Debug_counter > length(param.DEBUG_Count))
                        status = 17;
                        return;
                    else
                        Debug_counter = Debug_counter + 1;
                    end
                end
                i = rows(count);
                j = cols(count);
                y_current = squeeze(I{1}(i,j,:));
                w_current = squeeze(W{1}(i,j,:));
                x_indices = rows-i;
                y_indices = cols-j;
                d = (x_indices.^2+y_indices.^2);
                [~,idx] = sort(d);
                nearest_x = rows(idx(2:param.N_neighbors+1));
                nearest_y = cols(idx(2:param.N_neighbors+1));
                linear_indices = sub2ind([Nx,Ny], nearest_x, nearest_y);
                I_temp = reshape(I{1},[Nx*Ny,Nt]);
                W_temp = reshape(W{1},[Nx*Ny,Nt]);
                y_neighbors = I_temp(linear_indices, :);
                w_neighbors = W_temp(linear_indices, :);
                nearest_x_tshift = rows(idx(1:min(N_pixels,20)));
                nearest_y_tshift = cols(idx(1:min(N_pixels,20)));
                linear_indices_tshift = sub2ind([size(I{1}, 1), size(I{1}, 2)], nearest_x_tshift, nearest_y_tshift);
                y_avg_tshift = mean(I_temp(linear_indices_tshift, :));
                if (param.AmpEst_flag == 1)
                    Amp_temp = param.AmpMap(i,j);
                    Prior_Current{2} = [{1-Amp_temp},{Amp_temp}];
                end
                if ~isempty(param.PriorData)
                    Prior_temp = cell(param.order,1);
                    Prior_temp{2} = param.PriorData(i,j);
                    Prior_Current{1} = Prior_temp;
                end
                if (param.NeighborsNoise_flag == 1 || param.Method == 7)
                    NoiseVector = squeeze(Noise_calib(i,j,:));
                else
                    NoiseVector = [];
                end
                if (param.SweepCostfn == 1)
                    mydeconv_SweepCostfn(y_current, t_sig, IRF_current, t_IRF,param,NoiseVector);
                end
                [y_fit,Noise_var,amp,tau,chi,err_status,chi_vec,Ar_out,DC_out,tshift_out,err_str,param] = mydeconv8(y_current, w_current, t_sig{1}, IRF_current, t_IRF,param,NoiseVector,y_neighbors,w_neighbors,y_avg_tshift,Prior_Current);
                

                if (err_status == 0)
                    I_fit{count} = y_fit;
                    amp_map_vec{count} = amp;
                    Noise_Var_vec{count} = double(Noise_var);
                    tau_map_vec{count} = tau;
                    chi_map_vec{count}  = chi;
                    chi_vec_arr{count}  = chi_vec;
                    Ar_map_vec{count}  = Ar_out;
                    DC_map_vec{count}  = DC_out;
                    tshift_map_vec{count}  = tshift_out;
                else     %%% If we caught an error
                    I_fit{count} = [];
                    amp_map_vec{count}  = NaN;
                    Noise_Var_vec{count} = [];
                    tau_map_vec{count}  = NaN;
                    chi_map_vec{count}  = NaN;
                    Ar_map_vec{count}   = NaN;
                    DC_map_vec{count}  = NaN;
                    tshift_map_vec{count}  = NaN;
                    chi_vec_arr{count}  = [];
                end
                err_status_vec{count} = err_status;
            end
        end
        
        %%% Fill the output matrices
        try
            for count = 1:N_pixels
                i = rows(count);
                j = cols(count);
                tau_map(i,j,:) = tau_map_vec{count};
                if (param.Method ~= 20)
                    chi_map(i,j)  = chi_map_vec{count};
                    amp_map(i,j,:) = amp_map_vec{count};
                    Ar_map(i,j,:) = Ar_map_vec{count};
                    DC_map(i,j,:) = DC_map_vec{count};
                    tshift_map(i,j,:) = tshift_map_vec{count};
                end
            end
            Returned_data.tau_map = tau_map;
            Returned_data.amp_map = amp_map;
            Returned_data.chi_map = chi_map;  
        catch ME
            if (param.DEBUGGING == 1)
                 keyboard;
            else
                disp(strcat('Error in copying matrix! At count: ',num2str(count)))
                disp(ME.message)
                status =  7;
                return
            end
        end
        t_final = toc/60;
        fprintf(' (%.2f mins)',round(t_final,2))
    end
    if (param.Analysis_type == 3)

        %%%%% Common mask
        [Nx,Ny,~] = size(I{1});
        Mask_common = ones(Nx,Ny);
        Nt_all = 0;
        N_folders = size(t_sig,2);
        Nt = cell(N_folders,1);
        for idx_folder = 1:N_folders
            Nt{idx_folder} = size(I{idx_folder},3);
            Mask_common = Mask_common.*param_folder{idx_folder}.Mask;
            Nt_all = Nt_all+Nt{idx_folder};
        end
        %%%%% Extract Objects
        N_pixels_all = sum(Label_Mask(:));
        if (N_pixels_all == 0)
            Returned_data.flag_empty   = 1;
            return;
        end
        %%%%% Form noise matrix
        x_indices = (1:Nx).';
        y_indices = 1:Ny;
        W = cell(N_folders,1);
        for idx_folder = 1:N_folders
            %%% If spatial neighbor distance is needed
            [rows, cols] = find(Label_Mask > 0);
            N_pixels = length(rows);
            W_local = zeros(N_pixels, Nt{idx_folder});      % worker-safe sliced variable
            parfor idx_pixel = 1:N_pixels
                i = rows(idx_pixel);j = cols(idx_pixel);
                d = ((x_indices-i).^2+(y_indices-j).^2);
                [~,idx] = sort(d(:),'ascend');
                nn_lin = idx(1:param.N_neighbors);
                [i_near,j_near] = ind2sub([Nx,Ny],nn_lin);
                y_neighbors = I{idx_folder}(i_near,j_near,:);
                y_neighbors_avg = squeeze(mean(y_neighbors,[1,2]));
                y_neighbors_smoothed = smooth(y_neighbors_avg, 7);
                NoiseVector_empirical  = param.MCPNoise(cast(abs(y_neighbors_smoothed),"int16"));
                W_local(idx_pixel, :) = (1 ./ NoiseVector_empirical).';
            end
            W_temp = zeros(Nx,Ny,Nt{idx_folder});
            for p = 1:N_pixels
                W_temp(rows(p), cols(p), :) = reshape(W_local(p, :), 1, 1, Nt{idx_folder});
            end
            W{idx_folder} = W_temp;
        end

        if (param.cell_level_analysis == 1)
            %%%%% Concatenate in time.
            N_Object = length(unique(Label_Mask))-1;
            Returned_data = cell_level_analysis_v1(I,W,IRF,t_sig,t_IRF,param_folder,Label_Mask,N_Object);
            Returned_data.Type = "Exp";
            Returned_data.flag_empty   = 0;    
        elseif (param.ControlsCalib == 1)
            Returned_data = ControlsCalib(I,W,IRF,t_sig,t_IRF,param_folder);
        end
    end
    function updateProgress(~)
    persistent progress
    persistent last_print
    if isempty(progress)
        progress = 0;
        last_print = 0;
    end
    progress = progress + 1;
    percent_done = 100 * progress / N_pixels;
    if percent_done >= last_print + 0.1 || percent_done == 100
        if (last_print ~= 0)
            if (percent_done < 10)
                fprintf('\b\b\b\b');
            else
                fprintf('\b\b\b\b\b');     
            end
        end
        fprintf('%.1f%%', percent_done);
        last_print = percent_done;
    end
    if (percent_done == 100)
        progress = 0;
        last_print = 0;
    end
end
end

%% Extraction functions

function [freq0,delta_t,tmin,tmax,Gate,BINNING,MCP,SubDrk] = GetExperimentalParameters(parametersFile)
    A = fileread(parametersFile);
    A = A(~isspace(A)); % remove all spaces
    % need to extract experimental parameters
    str1 = 'TrigFreq:'; % for laser frequency
    str2 = 'Res:';  % for time steps
    str3 = 'On';    % for binning 
    str3_2 = 'Off';    % for binning 
    str4 = 'LaserTriggerRate='; % for laser frequency (if initial identifier fails)
    str5 = 'ScanRange:'; % for laser frequency (if initial identifier fails)
    str6 = 'CombHi';
    str7 = 'MCPGain:';
    str8 = 'Cam1:';
    % extract laser frequency
    indexFreq_1 = strfind(A,str1);
    indexFreq_2 = indexFreq_1+length(str1);
    indexFreq_3 = strfind(A,str4);
    indexFreq_4 = indexFreq_3+length(str4);
    
    
    freq0 = str2double(extractBefore(A(indexFreq_2:length(A)),'MHz'))*1E6;  % convert [MHz] to [Hz]
    if(isnan(freq0))
        % if the initial extraction fails, use the other identifier
        freq0 = str2double(extractBefore(A(indexFreq_4:length(A)),'MHz'))*1E6;  % convert [MHz] to [Hz]
    end
    % extract time step
    indexTime_1 = strfind(A,str2);
    indexTime_2 = indexTime_1+length(str2);
    indexTime_3 = strfind(A,str5)+length(str5);

    
    delta_t = str2double(extractBefore(A(indexTime_2:length(A)),'ps'))*1E-3;   % convert [ps] to [ns]
    tmin = str2double(extractBefore(A(indexTime_3:length(A)),'ps'))*1E-3;   % convert [ps] to [ns]\
    indexTime_4 = indexTime_3 + length(extractBefore(A(indexTime_3:length(A)),'ps'))+3;
    tmax = str2double(extractBefore(A(indexTime_4:length(A)),'ps'))*1E-3;   % convert [ps] to [ns]

    indexGate = strfind(A,str6)+length(str6);
    Gate = str2double(extractBefore(A(indexGate:indexGate+5),'ps'));

    indexGate = strfind(A,str7)+length(str7);
    MCP = str2double(extractBefore(A(indexGate:indexGate+5),'V'));

    indexSubDrk = strfind(A,str8);
    str_temp = A(indexSubDrk:indexSubDrk+20);
    indexSubDrk_2 = strfind(str_temp,'ms');
    str_temp2 = str_temp(indexSubDrk_2+2:indexSubDrk_2+3);
    if strcmp(str_temp2,'On')
        SubDrk = 1;
        indexBin_1 = strfind(A,str3);
        indexBin_2 = indexBin_1+length(str3);
        BINNING = str2double(extractBefore(A(indexBin_2:length(A)),'x'));
    else
        SubDrk = 0;
        indexBin_1 = strfind(A,str3_2);
        indexBin_2 = indexBin_1+length(str3_2);
        BINNING = str2double(extractBefore(A(indexBin_2:length(A)),'x'));
    end
end

function Noise_calib = Get_Noise_Neighbors(param,I)
    fprintf("Starting noise calibration... ");
    x_coord = 1:size(param.Mask,1);
    y_coord = 1:size(param.Mask,2);
    x_c_mat = repmat(reshape(x_coord,[length(x_coord),1]),1,length(y_coord));
    y_c_mat = repmat(y_coord,length(x_coord),1);
    N_t = length(t_sig);
    Noise_calib = zeros(size(I,1),size(I,2),N_t);
    fit_pix_x = [];
    fit_pix_y = [];
    i_image = sum(I,3);
    for q = 1:N_Obj
        %box each object
        x_c_mat_f = (Label_Mask == q).*x_c_mat;
        y_c_mat_f = (Label_Mask == q).*y_c_mat;
        x_c_image = nonzeros(x_c_mat_f);
        y_c_image = nonzeros(y_c_mat_f);
        xrange = min(x_c_image):max(x_c_image);
        yrange = min(y_c_image):max(y_c_image);
        %%% display box
        % figure(999 + q);
        % imagesc(i_image(xrange,yrange));
        %%% add included pixels to list for lifetime fitting
        fit_pix_x = [fit_pix_x;x_c_image];
        fit_pix_y = [fit_pix_y;y_c_image];

        %perform noise calculation and segmentation for box
        %segment intensity image w/ kmeans clustering
        image_seg = imsegkmeans(uint16(i_image(xrange,yrange)),param.KMeansSeg); 
        %%todo: make number of object tunable paramter
        %figure(100 + q);
        %imagesc(image_seg);

        x_idx = reshape(xrange,[length(xrange),1]);
        x_mat = repmat(x_idx,1,length(yrange));
        y_idx = yrange;
        y_mat = repmat(y_idx,length(xrange),1);

        parfor i = xrange
            for j = yrange
                region = image_seg(i-xrange(1)+1,j-yrange(1)+1);
                mask = (image_seg == region);
                x_indices = x_mat-i;
                y_indices = y_mat-j;
                distance = (x_indices.^2+y_indices.^2).*mask + (mask == 0).*999999999;
                [testmindist,idx] = sort(reshape(distance,[length(xrange)*length(yrange),1]));
                idx_x = mod(idx(1:25),length(xrange)) + length(xrange)*(mod(idx(1:25),length(xrange))==0);
                idx_y = ceil(idx(1:25)/length(xrange) + 0.0001);
                for p = 1:N_t
                    noise_points = zeros(1,25);
                    for r = 1:25
                        noise_points(1,r) = I(idx_x(r)+xrange(1)-1,idx_y(r)+yrange(1)-1,p);
                    end
                    %noise_calib(i,j,p) = var(noise_points);
                    stdnp = std(noise_points);
                    meannp = mean(noise_points);
                    %meannp = median(noise_points);
                    a2 = [];
                    for r = 1:25
                        %todo - add if statement back if working with
                        %real data
                        if abs(noise_points(1,r) - meannp) <= 3*(stdnp)
                            a2 = [a2;noise_points(1,r)];
                        end
                    end
                    if ~(isnan(var(a2)))
                        Noise_calib(i,j,p) = var(a2);
                    else
                        Noise_calib(i,j,p) = 0;
                    end
                end
            end
        end
    end
    % if (param.DEBUGGING == 1 && param.PLOT_FLAG == 1)
    %     noise_calib_img = sum(Noise_calib,3)/length(Noise_calib);
    %     figure(110);
    %     imagesc(noise_calib_img);
    % end
    fprintf("Done! ");
end


%% Saving/Plotting data
function Return_vector = plotObjLabels11_Labels(param,param_plot,path,Return_vector)
    %%% Set filenames and load results
    Excel_path = Get_savefile_name(param,path,"Excel");
    Results_filepath  = Get_savefile_name(param,path,"Results");
    Mask_filepath = Get_savefile_name(param,path,"Mask");
    Label_filepath = Get_savefile_name(param,path,"Label");
    fig_name_template = Get_savefile_name(param,path,"Template");
    flag_No_Mask = 0;
    Stain_flag = 0;
    if (~isfile(Mask_filepath))
        disp("Can't find Threshold Mask")
        flag_No_Mask = 1;
    else
        M = load(Mask_filepath);
        if (isfield(M,'Skip') && M.Skip == 1) || (isfield(M,'Skip_once') && M.Skip_once == 1)
            fprintf(" - Skipped!")
            Return_vector = struct("status",-1);
            return;
        else
            Mask = M.Mask;
        end
    end
    
    if(isfile(Label_filepath))
        load(Label_filepath,'Label_Mask','Type_Mask','N_Obj','-mat');
    end
    
    if(~isfile(Results_filepath))
        disp("Can't find Fitting Results: "+ string(Results_filepath))
        Return_vector = struct("status",0);
        return
    else
        if ismember(param.Analysis_type,[1,2])
            try
                load(Results_filepath,'Returned_data')
                I_avg = Returned_data.I_avg;
                Peak_Img = Returned_data.Peak_Img;
                tau_map = Returned_data.tau_map;
                amp_map = Returned_data.amp_map;
                chi_map = Returned_data.chi_map;
                HorizontalPixels = Returned_data.HorizontalPixels;
                VerticalPixels = Returned_data.VerticalPixels;
                time = Returned_data.time;
                if isfield(Returned_data,"Mask_used")
                    Mask = Returned_data.Mask_used;
                else
                    Mask = (tau_map > 0);
                end
                Mask = logical(Mask);
                if isfield(Returned_data,"Photon_Budget")
                    Photon_Budget = Returned_data.Photon_Budget;
                end
            catch ME
                err1 = ME.message;
                try
                    load(Results_filepath,'I_avg','Peak_Img','idx_vector','chi_seg','tau_map','amp_map','chi_map','HorizontalPixels','VerticalPixels','I_fit_seg','I_seg_actual','tau_seg','tau_dist_seg','time','Obj_Mask','chi_vector_map','PriorData')
                catch ME
                    disp(err1)
                    disp(ME.message)
                end
                % disp("Results loading error! (File actually exists) Check if Analysis_type is correct, is it really 1 or 2?")
            end
            if (param.Stain)
                if(isfield(Returned_data, 'StainImg'))
                    I_stain = Returned_data.StainImg;
                    fprintf(". Saved stain loaded.")
                    if (~isempty(I_stain))
                        Stain_flag = 1;
                    end            
                else
                    fprintf(". Stain not saved.")
                    if (param.current_stain_exists_flag)
                        try
                            N_files = length(param.current_stain_filenames);
                            for idx_stain = 1:N_files
                                Image_handle = Tiff(fullfile(fileparts(fileparts(path)),param.current_stain_filenames{idx_stain}),'r');  % read given frame
                                I_stain(:,:,idx_stain) = read(Image_handle);
                            end
                            fprintf(". Stain loaded.")
                            if (~isempty(I_stain))
                                Stain_flag = 1;
                            end
                        catch
                            fprintf(". Stain not loaded.")
                            Stain_flag = 0;
                        end
                    else
                        fprintf(". No stain found.")
                        Stain_flag = 0;
                    end
                end
                if (Stain_flag)
                    if (~isempty(I_stain))
                        I_stain = mean(I_stain,3);
                    end
                end
            end

        elseif ismember(param.Analysis_type,3)
            load(Results_filepath,'Returned_data')
            if (Returned_data.flag_empty   == 0)
                I_avg = Returned_data.I_avg;
                Peak_Img_arr = Returned_data.Peak_Img;
                Peak_Img = Peak_Img_arr{1};
                theta_map = Returned_data.theta_map;
                theta_avg = Returned_data.theta_avg;
                z_map = Returned_data.z_map;
                z_avg = Returned_data.z_avg;
                chi_map = Returned_data.chi_map;
                HorizontalPixels = Returned_data.HorizontalPixels;
                VerticalPixels = Returned_data.VerticalPixels;
                time = Returned_data.time;
                Obj_Mask = Returned_data.Obj_Mask;
                N_Obj = size(theta_avg,1);
                for idx_obj = 1:N_Obj
                    tau_map{idx_obj} = theta_map{1}{idx_obj}(2,:);
                end
                if(isfield(Returned_data, 'StainImg'))
                    %%% In Cell level analysis we have to have the
                    %%% staiin accounted for in the analysis for each
                    %%% cell earlier
                    I_stain = Returned_data.StainImg;
                    fprintf(". Saved stain loaded.")
                    if (~isempty(I_stain))
                        Stain_flag = 1;
                        stats = regionprops(Obj_Mask, Returned_data.Size_mask, 'MeanIntensity');
                        size_arr = [stats.MeanIntensity]';
                        N_Object_auto = N_Obj;
                    end    
                end
            else
                disp("Data skipped due to low number of pixels")
                return;
            end
        end
        try
            load(Results_filepath,'Noise_Var_vec');
        end
        % load(Results_filepath, '-regexp', '^(?!param)\w');
    end

    %%%% Edit the mask to remove saturated pixels..
    if iscell(Peak_Img)
        Mask(find(Peak_Img{1} == 4095)) = false;
        Peak_Img = Peak_Img{1};
    else
        Mask(find(Peak_Img == 4095)) = false;
    end
    N_pixels = sum(Mask(:));
    if (N_pixels == 0)
        return;
    end
    if ismember(param.Analysis_type,[1,2])
        if (param.Stain)
            if (Stain_flag) %%% Succeeded in loading it
                [Obj_Mask,celltype_array,cellsize_arr,Mask] = getHealthTypeAndSize(I_stain,Mask);
                N_Object_auto = length(cellsize_arr);
                % figure(1);imagesc(Mask)
                % figure(2);imagesc(Stain_mask)
                % figure(3);imagesc(Mask.*Stain_mask)
            end
        end
    end

    %%% Generating label mask
    if ismember(param.Analysis_type,[1,2])
        % Salem: Temporarily removed  for Tx data
        % Mask_temp  = bwareaopen(Mask,50,8);
        Mask_temp  = Mask;
        if (~param.Stain)
            if (~Stain_flag) %%% Succeeded in loading it
                [Obj_Mask,N_Object_auto] = bwlabel(Mask_temp);
            end
        end
    end

    %% Some Pre-calculations
    if ismember(param.Analysis_type,[1,2])
    
        Scalebar_length = param_plot.scale_bar;
        PixelSize = 11.22/param_plot.Magnification;
        Size_pixels = round(Scalebar_length/PixelSize);    
        tau_plot = tau_map;
        tau_plot(tau_plot == 0) = NaN;
        if iscell(Peak_Img)
            Peak_I_arr = Peak_Img{1}(Mask);
        else
            Peak_I_arr = Peak_Img(Mask);
        end
        dE = diff(param_plot.hbins_tau(1:2));
        hbins_chi = param_plot.hbins_chi;
        chi_avg = mean(chi_map(chi_map>0));
        chi_map (chi_map == 0) = NaN;
        [V_chi,E_chi] = histcounts(chi_map,param_plot.hbins_chi);
        dE_chi = diff(E_chi(1:2));
        W_chi = V_chi/sum(V_chi);
        chi_mean = sum(E_chi(1:end-1).*W_chi);
        chi_std = sqrt(sum(W_chi.*(E_chi(1:end-1) + dE_chi/2 - chi_mean).^2));
        
        Return_vector.chi_mean = chi_mean;
        Return_vector.chi_std = chi_std;
        Return_vector.W_chi = W_chi;
        Return_vector.hbins_chi = E_chi(1:end-1);
        
        E_temp = param_plot.hbins_tau;
        tau_bins = E_temp(1:end-1);
    
        g_mu = zeros(param.order,1);
        g_std = zeros(param.order,1);
        gw_mu = zeros(param.order,1);
        gw_std = zeros(param.order,1);
        tau_arr_all = cell(param.order,1);
        Peak_arr_all = cell(param.order,1);
        if (param.order > 1)
            idx_bad = [];
            tau_temp_all = [];
            amp_temp_all = [];
            idx_bad_all = [];
            for i = 1:param.order
                tau_temp = tau_plot(:,:,i);
                tau_temp  = tau_temp(:);
                amp_temp = amp_map(:,:,i);
                amp_temp = amp_temp(:);
                idx_bad = find(tau_temp<=param.lb+0.05 | tau_temp>=param.ub-0.1);
                if (param.Method == 11)
                    idx_bad = find(tau_temp<=0.52 | tau_temp>=5.5); 
                end
                idx_bad_all = [idx_bad_all;idx_bad];
                tau_nan_idx = isnan(tau_temp(:));
                amp_temp(tau_nan_idx) = NaN;
                tau_idx = ~isnan(tau_temp(:));
                tau_temp = tau_temp(tau_idx);
                amp_temp = amp_temp(tau_idx);   
                tau_temp_all = [tau_temp_all;tau_temp];
                amp_temp_all = [amp_temp_all;amp_temp];
                [V_temp,E_temp] = histcounts(tau_temp,param_plot.hbins_tau);             
                [V_amp_temp,E_amp_temp] = histcounts(amp_temp,param_plot.hbins_amp);
                [V_temp2, ~] = histwv(tau_temp, amp_temp, 0, 10, length(param_plot.hbins_tau));
                tau_dist{i} = V_temp/sum(V_temp);
                tau_mean_temp = sum(tau_bins.*tau_dist{i});
                tau_std_temp  = sqrt(sum(tau_dist{i}.*(tau_bins + dE/2 - tau_mean_temp).^2));
                tau_dist_weight_temp= V_temp2'/sum(V_temp2);
                tau_dist_weight{i}= tau_dist_weight_temp(1:end-1);
                tau_mean_temp_w = sum(tau_bins.*tau_dist_weight{i});
                tau_std_temp_w  = sqrt(sum(tau_dist_weight{i}.*(tau_bins + dE/2 - tau_mean_temp_w).^2));
        
                g_mu(i) = tau_mean_temp;
                g_std(i) = tau_std_temp;
                gw_mu(i) = tau_mean_temp_w;
                gw_std(i) = tau_std_temp_w;
    
                tau_arr = tau_map(Mask);
                [Peak_arr_all{i},idx_sorted] =  sort(Peak_I_arr);
                tau_arr_all{i} = tau_arr(idx_sorted);
            end
            N_bad_pixels = length(unique(idx_bad_all));
            tau_vector = param_plot.hbins_tau;
            Return_vector.SR = N_bad_pixels/sum(Mask(:));
        elseif (param.order == 1)
            tau_temp = tau_plot(:);
            tau_idx = ~isnan(tau_temp);
            tau_temp = tau_temp(tau_idx);    
            [V_temp,E_temp] = histcounts(tau_temp,param_plot.hbins_tau);             
            tau_dist = V_temp/sum(V_temp);
            tau_mean_temp = sum(tau_bins.*tau_dist);
            tau_std_temp  = sqrt(sum(tau_dist.*(tau_bins + dE/2 - tau_mean_temp).^2));
            g_mu = tau_mean_temp;
            g_std = tau_std_temp;
            IC = [g_mu,g_std];
            [g_amp,gw_mu,gw_std,fit_hist,x] = FitHistogram_tau(tau_bins,tau_dist,param,IC);
    
            tau_arr = tau_map(Mask);
            chi_arr = chi_map(Mask);
            [Peak_arr_all,idx_sorted] =  sort(Peak_I_arr);
            tau_arr_all = tau_arr(idx_sorted);
            Chi_arr_all = chi_arr(idx_sorted);
        end
        Return_vector.MCP = param.MCP;
        Return_vector.tau = g_mu;
        Return_vector.std = g_std;
        Return_vector.tauw = gw_mu;
        Return_vector.stdw = gw_std;    
        Return_vector.tau_arr = tau_arr_all;
        Return_vector.Peak_arr = Peak_arr_all;
        try
            Return_vector.Chi_arr = Chi_arr_all;
        catch
        end
        Return_vector.tau_map = tau_map;
        Return_vector.Peak_Img = Peak_Img;
        Return_vector.Mask = Mask;
        Return_vector.path = path;
        Return_vector.time = time;
        if exist('Noise_Var_vec','var')
            Return_vector.Noise_var_vec = Noise_Var_vec;
        else
            Return_vector.Noise_var_vec = [];
        end
    elseif ismember(param.Analysis_type,3)
        Return_vector = Returned_data;
    end
    %% Set dimensions based on Objects Number
    %%%% This is skipped because the Label file idea is not used at all and
    %%%% will potentially be removed or modified later
    % if (isfile(Label_filepath))
    %     N_Obj_actual  =  0;
    %     for idx_Obj = 1:N_Obj   %%% Quick loop to count actual labels used
    %         FullMask_temp = Label_Mask;
    %         FullMask_temp ((Label_Mask ~= idx_Obj) & (Label_Mask ~= 0)) = 1;
    %         FullMask_temp (Label_Mask == idx_Obj) = 2;
    %         N_Obj_actual = N_Obj_actual + 1;
    %         if (sum(FullMask_temp == 2) == 0)
    %               N_Obj_actual = N_Obj_actual - 1;
    %         end
    %     end
    % 
    %     if (N_Obj_actual == 2)
    %         fig1_Pos = [50 50 1100  500];
    %         ax1_pos  = [-0.25 0 1 1];
    %         ax_hist  = @(i) [0.6 0.6-(i-1)*0.5 0.2 0.25];
    %         leg_pos = [0.85 0.5 0.28 0.28];
    %     elseif (N_Obj_actual <= 3)
    %         fig1_Pos = [50 50 1100 800];
    %         ax1_pos  = [0 0.25 0.6 0.6];
    %         ax_hist  = @(i) [0.62 0.72-(i-1)*0.33 0.2 0.2];
    %         leg_pos = [0.9 0.4 0.28 0.28];
    %     elseif (N_Obj_actual > 3 && N_Obj_actual <= 6)
    %         fig1_Pos = [50 50 1450 800];
    %         ax1_pos  = [-0.1 0.25 0.6  0.6];
    %         ax_hist = @(i) [0.42+0.25*floor(mod((i-1)/3,3)) 0.72-(i-1-3*floor(mod((i-1)/3,3)))*0.33 0.2 0.2];
    %         leg_pos = [0.9 0.4 0.28 0.28];
    %     elseif (N_Obj_actual > 6 && N_Obj_actual <= 9)
    %         fig1_Pos = [50 50 1800 800];
    %         ax1_pos  = [-0.15 0.25 0.6  0.6];
    %         ax_hist  = @(i) [0.34+0.22*floor(mod((i-1)/3,3)) 0.72-(i-1-3*floor(mod((i-1)/3,3)))*0.33 0.15 0.2];
    %         leg_pos  = [0.1 0 0.28 0.28];
    %     end
    % end
    
    %% Plot the main figure showing the Objects in separate colors
    if (param_plot.plot_ObjHist == 1)
        if (isfile(Label_filepath))
            gcf = figure(1);
            Mycolormap  = [[0 0 1];[1 0 0];[0 1 0];[1 0 1];[1 1 0];[0 1 1];[1 1 1];[0.8 0.8 1];[0.5 1 0.5]]; %% Make this bigger
            gcf.Position = fig1_Pos;
            ax1 = axes(gcf,'Position',ax1_pos);
            ax2 = copyobj(ax1,gcf);
            ax3 = copyobj(ax1,gcf);
            imagesc(ax1,Peak_Img/max(max(Peak_Img))); %% Plot the image
            clim(ax1,[param_plot.clim_I_min,param_plot.clim_I_max])
            xmin = 100;
            xmax = 800;
            xlim([xmin ,xmax])
            
            imagesc(ax2,Label_Mask, 'AlphaData', 0.5*(Label_Mask > 0));
            if (N_Obj_actual > 1)
                clim(ax2,[1,N_Obj_actual])
            end
            ax2.UserData = linkprop([ax1,ax2],{'Position','InnerPosition','DataAspectRatio','xtick','ytick','ydir','xdir','xlim','ylim'});
            ax3.UserData = linkprop([ax1,ax3],{'Position','InnerPosition','DataAspectRatio','xtick','ytick','ydir','xdir','xlim','ylim'});
            
            title('Objects maps','Fontsize',20);    
            colormap(ax1,'gray')
            colormap(ax2,Mycolormap(1:N_Obj_actual,:))
            
            if (param_plot.scale_Left_or_right == 1)
                x_vector = [xmin+20 xmin+Size_pixels+20];   % for plotting
            else
                x_vector = [xmax-20-Size_pixels xmax-20];   % for plotting
            end
            y_vector = [size(Label_Mask,1)-20 size(Label_Mask,1)-20];   % for plotting
            line(ax3,x_vector,y_vector,'LineWidth',param_plot.scale_bar_Linewidth,'Color',[0.99 0.99 0.99]);  % plot line
            set(ax1,'XColor', 'none','YColor','none')
            set(ax2,'XColor', 'none','YColor','none')
            set(ax3,'XColor', 'none','YColor','none')
            ax2.Visible = 'off';
            ax3.Visible = 'off';
            ax = cell(N_Obj_actual,1);
            i_actual_Obj = 0;
            for idx_Obj = 1:N_Obj   %%  for each label separately
                FullMask_temp = Label_Mask;
                FullMask_temp ((Label_Mask ~= idx_Obj) & (Label_Mask ~= 0)) = 1;
                FullMask_temp (Label_Mask == idx_Obj) = 2;
                i_actual_Obj = i_actual_Obj + 1;
                if (sum(FullMask_temp(:) == 2) == 0)
                      disp(strcat("Not enough points for label: ",num2str(idx_Obj)));
                      i_actual_Obj = i_actual_Obj - 1;
                      continue;
                end
                %%%  Objects Histogram
                ax{idx_Obj} = axes('Position',ax_hist(i_actual_Obj));hold on;
                box(ax{idx_Obj},'on')
                set(ax{idx_Obj},'LineWidth',2,'XColor',Mycolormap(idx_Obj,:),'YColor',Mycolormap(idx_Obj,:));    
                Mask_Temp = Label_Mask;
                Mask_Temp(Mask_Temp ~= idx_Obj) = 0;
                Mask_Temp(Mask_Temp == idx_Obj) = 1;
                tau_Obj = tau_plot;
                tau_Obj(~Mask_Temp) = NaN;
                [V,E] = histcounts(tau_Obj,param_plot.hbins_tau);
                W = V/sum(V);
                dE = diff(E(1:2));
                tau_mean = sum((E(1:end-1)+dE/2).*W);
                tau_std = sqrt(sum(W.*(E(1:end-1) + dE/2 - tau_mean).^2));
            
                chi_Obj = chi_map;
                chi_Obj(~Mask_Temp) = NaN;    
                [V_chi,E_chi] = histcounts(chi_Obj ,param_plot.hbins_chi);
                W_chi = V_chi/sum(V_chi); 
                dE_chi = diff(E_chi(1:2));
                chi_mean = sum((E_chi(1:end-1)+dE_chi/2).*W_chi);
                chi_std = sqrt(sum(W_chi.*(E_chi(1:end-1) + dE_chi/2 - chi_mean).^2));
            
                plot(ax{idx_Obj},E(1:end-1),smooth(W),'color',Mycolormap(idx_Obj,:),'LineWidth',2);
                set(ax{idx_Obj},'XLim',[param_plot.hist_xlim_min,param_plot.hist_xlim_max],'fontweight','bold')
                hx = xlabel('Lifetime (ns)','FontSize',15);
                hy = ylabel('Distribution','FontSize',15);
                title_str = strcat("Lifetime mean: ",num2str(tau_mean));
                title_str2 = strcat("\chi^2 mean: ",num2str(chi_mean));
                title(ax{idx_Obj},[title_str;title_str2],'fontsize',13,'color',Mycolormap(idx_Obj,:))         
            end
            
            ax_temp = axes('Position',leg_pos);hold on;
            ax_temp .Visible = 'off';
            set(ax_temp ,'XColor', 'none','YColor','none')
            
            fig_name = get_file_name('ObjHist',fig_name_template,'.fig');
            png_name = get_file_name('ObjHist',fig_name_template,param_plot.figformat); 
            saveas(figure(1),fig_name)
            saveas(figure(1),png_name)
            close(figure(1));
        else
            disp("Can't find Label Mask")
        end
    end
    
    %% Plotting all the information
    if (param_plot.plot_mix == 1)
        gcf = figure(1);
        clf(gcf);    
        if (param.order == 1)
            ax1 = axes(gcf,'InnerPosition',[0.1300 0.1100 0.72 0.8150] );
            subplot(1,4,2);hold on;
            histogram(tau_plot,param_plot.hbins_tau);
            xlabel('Lifetime (ns)','fontsize', 15,'fontweight', 'bold');
            ylabel('Count','fontsize', 15,'fontweight', 'bold');
            title(['Mean = ' num2str(round(tau_mean,2)) ', Std = ' num2str(round(tau_std,3))],'fontsize', 14,'fontweight', 'bold')
            xline(tau_mean,'Color','red','LineWidth',1.5)
            xlim([max(0,tau_mean-10*tau_std),tau_mean+10*tau_std]);
            pbaspect([1 1 1]);
            IC = [max(V),tau_mean,tau_std];
            if (param_plot.Fit_Histogram == 1)
                pdf_temp = normpdf(E,IC(2),IC(3));
                pdf_temp = IC(1)*pdf_temp/max(pdf_temp);
                plot(E+dE/2,pdf_temp,'r--','LineWidth',2)
            end
            ax2 = subplot(1,4,1);
            imagesc(ax2,HorizontalPixels, VerticalPixels, tau_plot);
            colormap(ax2,'turbo')
            colorbar;
            clim([max(0,tau_mean-3*tau_std),tau_mean+3*tau_std]);
            title('Lifetime (ns)','fontsize', 14,'fontweight', 'bold');
            pbaspect([1 1 1]);
            ax3 = subplot(1,4,3);
            imagesc(HorizontalPixels, VerticalPixels, chi_map);
            colormap(ax3,'turbo')
            title(['\chi^2_{r} (Mean: ' num2str(round(chi_mean,2)) ', Std:' num2str(round(chi_std,1)) ')'],'fontsize', 14,'fontweight', 'bold');
            clim([max(0,chi_mean-1*chi_std),chi_mean+1*chi_std])
            colorbar;
            pbaspect([1 1 1]);
            
            ax4 = subplot(1,4,4);
            histogram(chi_map(:),param_plot.hbins_chi);
            xlim([chi_mean-3*chi_std,chi_mean+3*chi_std])
            colormap(ax3,'turbo')
            title('\chi^2_{r} (Histogram)','fontsize', 14,'fontweight', 'bold');
            pbaspect([1 1 1]);
    
            set(gcf,'position',[10,100,1400,400])
    
            
    
        elseif (param.order > 1)
            colormap jet
            IC = [];
            
            for i = 1:param.order
                subplot(param.order,5,2+5*(i-1));hold on;
                % h = histogram([tau_plot(:,:,i)],param_plot.hbins_tau);
    
                tau_temp = tau_plot(:,:,i);
                amp_temp = amp_map(:,:,i);
                tau_temp = tau_temp(:);
                amp_temp = amp_temp(:);
                tau_idx = ~isnan(tau_temp(:));
                tau_temp = tau_temp(tau_idx);
                amp_temp = amp_temp(tau_idx);
                
                % [V_temp, ~] = histwv(tau_temp, amp_temp, 0, 10, length(param_plot.hbins_tau));
                [V_temp, ~] = histcounts(tau_temp,param_plot.hbins_tau);
    
                bar(param_plot.hbins_tau(1:end-1),V_temp,'k','linewidth',2)
                E_temp = param_plot.hbins_tau(1:end-1);
                W_temp = (V_temp/sum(V_temp));
                tau_mean_temp = sum(E_temp.*W_temp);
                tau_std_temp  = sqrt(sum(W_temp.*(E_temp + dE/2 - tau_mean_temp).^2));
                xlabel('Lifetime (ns)','fontsize', 13,'fontweight', 'bold');
                ylabel('Count');
                title(['Mean = ' num2str(round(tau_mean_temp,1)) ',Std: ' num2str(round(tau_std_temp,1))])
                xline(tau_mean_temp,'Color','red','LineWidth',1.5)
                xlim([max(0,tau_mean_temp-3*tau_std_temp),tau_mean_temp+3*tau_std_temp]);
                % ylim([0,1.2*max(V_temp)])
                pbaspect([1 1 1]);
                IC = [IC max(V_temp),tau_mean_temp,tau_std_temp];
                
                subplot(param.order,5,1+5*(i-1));
                imagesc(HorizontalPixels, VerticalPixels, tau_plot(:,:,i));
                colorbar;
                clim([max(0,tau_mean_temp-3*tau_std_temp),tau_mean_temp+3*tau_std_temp]);
                title(strcat('lifetime ',num2str(i),' (ns)'));
                set(gcf,'position',[100,100,1200,400*param.order])
                pbaspect([1 1 1]);
            end                      
            subplot(param.order,5,3);hold on    
            
            %%% Weighted histogram here    
            tau_temp = tau_plot(:);
            amp_temp = amp_map(:);
            tau_idx = ~isnan(tau_temp(:));
            tau_temp = tau_temp(tau_idx);
            amp_temp = amp_map(tau_idx);                
            [V_temp, ~] = histcounts(tau_temp,param_plot.hbins_tau);
            bar(param_plot.hbins_tau(1:end-1),V_temp,'k','linewidth',2)
            title('Merged Histograms')
            pbaspect([1 1 1]);
    
    
            subplot(param.order,5,8);hold on    
            
            [V, ~] = histwv(tau_temp, amp_temp, 0, 10, length(param_plot.hbins_tau));
            bar(param_plot.hbins_tau,V,'k','linewidth',2)
            title('Weighted Merged Histograms')
    
            if (param_plot.Fit_Histogram == 1)
                tau_vector = param_plot.hbins_tau;
                [g_amp,g_mu,g_std,fit_hist,x] = FitHistogram_tau(param_plot.hbins_tau,V,param,IC);
                plot(tau_vector+dE/2,sum(V)*fit_hist/sum(fit_hist),'r--','linewidth',2);
                str_list = cell(param.order,1);
                for i = 1:param.order
                    str_list{i} = ['Mean = ' num2str(round(g_mu(i),2)) 'ns , Std: ' num2str(round(g_std(i),2)),' ns'];
                end
            title(str_list)
            end
            xlabel('Lifetime (ns)','fontsize', 13,'fontweight', 'bold');
            ylabel('Count');
            xlim([max(0,tau_mean-3*tau_std),tau_mean+3*tau_std]);
            pbaspect([1 1 1]);
            subplot(param.order,5,4)
            imagesc(HorizontalPixels, VerticalPixels, chi_map);
            colorbar;
            clim([max(0,chi_mean-1*chi_std),chi_mean+1*chi_std])
            title(['\chi^2_{r} (Mean: ' num2str(round(chi_mean,1)) ', Std:' num2str(round(chi_std,1)) ')'],'fontsize', 15,'fontweight', 'bold');
            pbaspect([1 1 1]);
            set(gcf,'position',[100,0,1200,800])
    
            subplot(param.order,5,9);
            histogram(chi_map(:),param_plot.hbins_chi);
            xlim([chi_mean-3*chi_std,chi_mean+3*chi_std])
            title('\chi^2_{r} (Histogram)','fontsize', 14,'fontweight', 'bold');
            pbaspect([1 1 1]);
    
            % set(gcf,'position',[10,100,1400,400])
    
            subplot(param.order,5,10);
            Peak_Img_temp = Peak_Img;
            Peak_Img_temp(Mask == 0) = NaN;
            [Peak_Img_sorted,Peak_idx] = sort(Peak_Img_temp(:));
            chi_temp = chi_map(:);
            chi_sorted = chi_temp(Peak_idx);
            plot(Peak_Img_sorted,chi_sorted);
    
            title('\chi^2_{r} vs intensity','fontsize', 14,'fontweight', 'bold');
            pbaspect([1 1 1]);
    
            
        end
        fig_name = get_file_name('Mix',fig_name_template,'.fig');
        png_name = get_file_name('Mix',fig_name_template,param_plot.figformat); 
        saveas(figure(1),fig_name)
        saveas(figure(1),png_name)
        close(figure(1));
    end
    
    %% Plot Chi Map
    if (param_plot.plot_chi == 1)
        gcf = figure(1);
        clf(gcf);
        ax = axes(gcf);
        xmin = param_plot.Img_xlim_min;
        xmax = param_plot.Img_xlim_max;
        ymin = param_plot.Img_ylim_min;
        ymax = param_plot.Img_ylim_max;
        gcf.Position = [0,0,xmax-xmin,ymax-ymin];
        ax.InnerPosition = [0,0,0.9,0.9];
        imagesc(ax,HorizontalPixels, VerticalPixels, chi_map);
        colormap jet
        colorbar;
        clim(ax,[max(0,chi_mean-1*chi_std),chi_mean+1*chi_std])
        % title(ax,['\chi^2_{r} (Mean: ' num2str(round(chi_mean,1)) ', Std:' num2str(round(chi_std,1)) ')'],'fontsize', 15,'fontweight', 'bold');
        set(ax,'XColor', 'none','YColor','none')
        xlim([xmin,xmax])
        ylim([ymin,ymax])
        pbaspect(ax,[xmax-xmin ymax-ymin 1])
        h2 = colorbar(ax);
        h2.Title.String = '\chi^2';
        h2.Location = 'eastoutside';
        h2.FontSize = 14;
        fig_name = get_file_name('ChiMap',fig_name_template,'.fig');
        png_name = get_file_name('ChiMap',fig_name_template,param_plot.figformat);
        saveas(figure(1),png_name)
        saveas(figure(1),fig_name)
        close(figure(1));
        %%
    %%%%%%%%%%%%%%%%%%% 
        % gcf = figure(1);
        % clf(gcf);
        % hold on;
        % % ax = axes(gcf);
        % [V,E] = histcounts(chi_map(:),param_plot.hbins_chi);
        % W_chi = V/sum(V);
        % bar(E(1:end-1),W_chi)
        % 
        % [chi_mean,chi_std,exp_dist] = FitHistogram_chi(E(1:end-1),W_chi,chi_mean,chi_std);
        % plot(E(1:end-1),exp_dist,'k--','linewidth',2.5)
        % xlim([chi_mean-3*chi_std,chi_mean+3*chi_std])
        % xlabel("\chi^2",'fontsize',15,'fontweight','bold')
        % ylabel("Distribution",'fontsize',15,'fontweight','bold')
        % xline(chi_mean,'Color','red','LineWidth',1.5)
        % legend('\chi^2_{r}','Gaussian fit','fontsize', 12,'fontweight', 'bold');
        % set(gca,'FontSize',15)
        % fig_name = get_file_name('ChiDist',fig_name_template,'.fig');
        % png_name = get_file_name('ChiDist',fig_name_template,param_plot.figformat);
        % saveas(figure(1),png_name)
        % saveas(figure(1),fig_name)
        % close(figure(1));
    %%%%%%%%%%%%%%%%%%%%%
    %%
        % gcf = figure(1);
        % clf(gcf);
        % ax = axes(gcf);   
        % Peak_Img_temp = Peak_Img;
        % Peak_Img_temp(Mask == 0) = NaN;
        % [Peak_Img_sorted,Peak_idx] = sort(Peak_Img_temp(:));
        % chi_temp = chi_map(:);
        % chi_sorted = chi_temp(Peak_idx);
        % I_count_arr = 1:20:3000;
        % chi_mean_arr = zeros(size(I_count_arr));
        % for idx_count = 1:length(I_count_arr)
        %     idx_peak = find(abs(Peak_Img_sorted- I_count_arr(idx_count)) < 40);
        % 
        %     if ~isempty(idx_peak)
        %         chi_std_arr = chi_sorted(idx_peak,:);
        %         [V_std,E_std] = histcounts(chi_std_arr,param_plot.hbins_chi);
        %         W = V_std/sum(V_std);
        %         E_std = E_std(1:end-1);
        %         chi_mean_arr(idx_count) = sum(E_std.*W);      
        %     else
        %         chi_mean_arr(idx_count) = NaN;
        %     end
        % end
        % 
        % plot(Peak_Img_sorted,chi_sorted); hold on
        % yline(1,'Color','k','LineWidth',2);
        % plot(I_count_arr,smooth(chi_mean_arr,50),'Color','red','LineWidth',2)
        % % plot(Peak_Img_sorted,smooth(chi_sorted,10),'r','linewidth',2); hold on
        % xlim([min(Peak_Img_temp(:)),max(Peak_Img_temp(:))])
        % ylim([0,chi_mean+20*chi_std])
        % xlabel("Intensity (Counts)",'fontsize',15,'FontWeight','bold')
        % ylabel("\chi^2",'fontsize',15,'FontWeight','bold')
        % legend('\chi^2','Ideal Mean','Estimated Mean','fontsize',12,'fontweight','bold')
        % fig_name = get_file_name('ChiInt',fig_name_template,'.fig');
        % png_name = get_file_name('ChiInt',fig_name_template,param_plot.figformat);
        % saveas(figure(1),png_name)
        % saveas(figure(1),fig_name)
        % close(figure(1));
    
        %%
        % gcf = figure(1);
        % clf(gcf);
        % ax = axes(gcf);   
        % I_temp = reshape(I,[size(I,1)*size(I,2),size(I,3)]);
        % Mask_temp = reshape(Mask,[size(I,1)*size(I,2),1]);
        % I_temp = I_temp(Mask_temp,:);
        % chi_temp = reshape(chi_vector_map,[size(I,1)*size(I,2),size(I,3)]);
        % chi_temp = chi_temp(Mask_temp,:);
        % 
        % [I_sorted,Peak_idx] = sort(I_temp (:));
        % chi_sorted = chi_temp(Peak_idx);
        % I_count_arr = 1:20:3000;
        % chi_mean_arr = zeros(size(I_count_arr));
        % for idx_count = 1:length(I_count_arr)
        %     idx_peak = find(abs(I_sorted- I_count_arr(idx_count)) < 40);
        % 
        %     if ~isempty(idx_peak)
        %         chi_std_arr = chi_sorted(idx_peak,:);
        %         [V_std,E_std] = histcounts(chi_std_arr,param_plot.hbins_chi);
        %         W = V_std/sum(V_std);
        %         E_std = E_std(1:end-1);
        %         chi_mean_arr(idx_count) = sum(E_std.*W);      
        %     else
        %         chi_mean_arr(idx_count) = NaN;
        %     end
        % end
        % 
        % plot(I_sorted,chi_sorted); hold on
        % yline(1,'Color','k','LineWidth',2);
        % plot(I_count_arr,smooth(chi_mean_arr,50),'Color','red','LineWidth',2)
        % % plot(Peak_Img_sorted,smooth(chi_sorted,10),'r','linewidth',2); hold on
        % % xlim([min(Peak_Img_temp(:)),max(Peak_Img_temp(:))])
        % % ylim([0,chi_mean+20*chi_std])
        % % ylim([0,10])
        % xlabel("Intensity (Counts)",'fontsize',15,'FontWeight','bold')
        % ylabel("\chi^2",'fontsize',15,'FontWeight','bold')
        % legend('\chi^2','Ideal Mean','Estimated Mean','fontsize',12,'fontweight','bold')
        % fig_name = get_file_name('ChiInt',fig_name_template,'.fig');
        % png_name = get_file_name('ChiInt',fig_name_template,param_plot.figformat);
        % saveas(figure(1),png_name)
        % saveas(figure(1),fig_name)
        % close(figure(1));
    end
    
    %% Plotting Histogram 
    if (param_plot.plot_hist == 1)      %%% Histogram plot only
        if (param_plot.Stack_figure ~= 1 && param.PLOT_FLAG == 1)
            figure(param_plot.fig_idx);
            clf(gcf);
        end
        [colorstr,str_leg,str_linestyle,str_MarkerStyle] = Get_Method_info(param);
        if (param.order > 1)
            tau_dist_all = zeros(1,length(tau_bins));
            for i = 1:param.order
                figure(param_plot.fig_idx);
                % if (i == 1)
                %     colorstr = 'b';
                % elseif (i == 2)
                %     colorstr = 'r';
                % end
                hold on;
                % N_line = 5;
                tau_dist_all = tau_dist_all + tau_dist_weight{i};
                % tau_dist_all = tau_dist_all + tau_dist{i};
                % if strcmp(param_plot.plot_hist_style,'bar')
                %     % bar(tau_bins,tau_dist{i} ,'EdgeColor',colorstr,'linewidth',1,"DisplayName",strcat("\tau_{",num2str(i),"} Histogram"),'EdgeAlpha',0.5)
                %     bar(tau_bins,tau_dist{i},'linewidth',1,"DisplayName",strcat("\tau_{",num2str(i),"} Histogram"),'EdgeAlpha',0.5,'EdgeColor',colorstr,'FaceColor',colorstr)
                % elseif strcmp(param_plot.plot_hist_style,'plot')
                %     plot(tau_bins,smooth(tau_dist{i},10),'linewidth',3,"DisplayName",strcat("\tau_{",num2str(i),"} Histogram"),'Color',colorstr)
                % end
                % line_tau = -0.2:0.025:0.2;
                % line_tau = line_tau + param_plot.delta_line;
                % plot(tau_bins,tau_dist{i} ,'Color',colorstr,'linewidth',2,"DisplayName",strcat("\tau_{",num2str(i),"} Histogram"))
                % line([param.NoiseActual_params(i) param.NoiseActual_params(i) param.NoiseActual_params(i) param.NoiseActual_params(i) param.NoiseActual_params(i)],[0 0.0025 0.005 0.0075 0.01], 'Color', colorstr,'Marker','o', 'LineWidth', 2,"DisplayName",strcat("\tau_{",num2str(i),"} true value")); % Plot the vertical line
                % line(param.NoiseActual_params(i)*ones(N_line,1),linspace(0,0.5,N_line), 'Color', [0 0 0], 'LineWidth', 2,"DisplayName", strcat("\tau_{", num2str(i), "} true value"),'LineStyle','--'); % Plot the vertical line
                % line(gw_mu(i)*ones(length(line_tau),1),line_tau, 'Color', colorstr, 'LineWidth', 2,"DisplayName", strcat("\tau_{", num2str(i), "} Avg"),'LineStyle','--',"Marker",str_MarkerStyle,'MarkerSize',10); % Plot the vertical line
                % line([gw_mu(i)   gw_mu(i) gw_mu(i) gw_mu(i)  gw_mu(i)], 0.001+[0 0.0025 0.005 0.0075 0.01], 'Color', colorstr, 'LineStyle','--', 'LineWidth', 2,"DisplayName",strcat("\tau_{",num2str(i),"} Estimated Value"),'Marker','x'); % Plot the vertical line
                % line([gw_mu(i)   gw_mu(i) gw_mu(i) gw_mu(i)  gw_mu(i)], 0.001+[0 0.0025 0.005 0.0075 0.01], 'Color', 'k', 'LineStyle','--', 'LineWidth', 2,"DisplayName",strcat("\tau_{",num2str(i),"} Estimated Value")); % Plot the vertical line
                xlabel("\tau (ns)")
                ylabel("Distribution (a.u.)")
                xxx = gca;
                xxx.FontSize = 14;
                % 
                % figure(5);hold on;
                % if (i == 1)
                %     clf;
                % end
                % plot(E_amp_temp(1:end-1),V_amp_temp./sum(V_amp_temp),'linewidth',2,"DisplayName",strcat("Lifetime ",num2str(i)))
                % title("Amplitude")
            end
            if strcmp(param_plot.plot_hist_style,'bar')
                bar(tau_bins,tau_dist_all/sum(tau_dist_all),'linewidth',1,"DisplayName",strcat("\tau_{",num2str(i),"} Histogram"),'EdgeAlpha',0.5,'EdgeColor',colorstr,'FaceColor',colorstr)
            elseif strcmp(param_plot.plot_hist_style,'plot')
                plot(tau_bins,smooth(tau_dist_all/sum(tau_dist_all),10),'linewidth',3,"DisplayName",strcat("\tau_{",num2str(i),"} Histogram"),'Color',colorstr)
            end

    
            %%% Plotting the CRB distribution. needs to be manually input
            if (param_plot.CRB_on_LT == 1)
                tau1 = 1;
                tau2 = 2.2;
                CRB1 = GetCRB(g_mu,param);
                CRB2 = GetCRB(flip(g_mu),param);
        
                pdf1 = normpdf(x,g_mu(1),CRB1);
                pdf2 = normpdf(x,g_mu(2),CRB2);
                plot(x,g_amp(1)*pdf1/sum(pdf1),'c','LineWidth',2,"DisplayName","CRB dist")
                plot(x,g_amp(2)*pdf2/sum(pdf2),'m','LineWidth',2)    
                str_list = cell(param.order,1);
                for i = 1:param.order
                    str_list{i} = ['Mean = ' num2str(round(g_mu(i),2)) 'ns , Std: ' num2str(round(g_std(i),2)),' ns'];
                end
                title(str_list)
            end
        elseif (param.order == 1)
            figure(param_plot.fig_idx+17);hold on;
            plot(tau_bins,smooth(tau_dist,5),'linewidth',2,'Color',colorstr,"DisplayName","\tau Histogram",'LineStyle','-')
            % line(param.NoiseActual_params(1)*ones(1,param_plot.N_points_marker),linspace(0,1.2*max(W_temp),param_plot.N_points_marker), 'Marker', 'o', 'LineWidth', 2,"DisplayName", strcat("\tau_{", num2str(i), "} true value")); % Plot the vertical line
            % line(g_mu*ones(1,param_plot.N_points_marker), linspace(0,1.2*max(tau_dist),param_plot.N_points_marker),'Color',colorstr, 'LineStyle','--', 'LineWidth', 2,"DisplayName",strcat("\tau_{",num2str(i),"} Estimated Value"),'Marker','x'); % Plot the vertical line
            xlabel("Lifetime (ns)")
            ylabel("Distribution")
            xxx = gca;
            xxx.FontSize = 14;
        end
        grid off;
        yticklabels([])    
    
        %% Histogram of Ar map
    
        % [V_Ar,E_Ar] = histcounts(Ar_map(:),0.6:0.01:1.6);
        % figure(8);hold on;
        % plot(E_Ar(2:end),V_Ar,'k','linewidth',2);
        % 
        % 
        % tau_temp = tau_plot(:,:,1);
        % DC_temp = DC_map(:);
        % tau_idx = ~isnan(tau_temp(:));
        % tau_temp = tau_temp(tau_idx);
        % DC_temp = DC_temp(tau_idx);  
        % 
        % [V_DC,E_DC] = histcounts(DC_temp(:),0:0.5:500);
        % figure(9);hold on;
        % plot(E_DC(2:end),V_DC,'k','linewidth',2);
        % 
        % tau_temp = tau_plot(:,:,1);
        % tshift_temp = tshift_map(:);
        % tau_idx = ~isnan(tau_temp(:));
        % tau_temp = tau_temp(tau_idx);
        % tshift_temp = tshift_temp(tau_idx);  
        % [V_tshift,E_tshift] = histcounts(tshift_temp(:),-0.5:0.01:0.5);
        % figure(10);hold on;
        % plot(E_tshift(2:end),V_tshift,'k','linewidth',2);
    
    
        % keyboard
    
        % if (param_plot.Stack_figure ~= 1)
        %     fig_name = get_file_name('Hist',fig_name_template,'.fig');
        %     png_name = get_file_name('Hist',fig_name_template,param_plot.figformat);
        %     saveas(figure(param_plot.fig_idx),fig_name)
        %     saveas(figure(param_plot.fig_idx),png_name)
        %     close(figure(param_plot.fig_idx));
        % end
    end
    
    %% Save Lifetime Images separately
    if (param_plot.plot_PeakLT == 1)
        gcf = figure(4);
        ax1 = axes(gcf);hold on;
        xmin = param_plot.Img_xlim_min;
        xmax = param_plot.Img_xlim_max;
        ymin = param_plot.Img_ylim_min;
        ymax = param_plot.Img_ylim_max;
    
        gcf.Position = [0,0,xmax-xmin,ymax-ymin];
        ax1.InnerPosition = [0,0,1,1];
        ax2 = copyobj(ax1,gcf);
        ax3 = copyobj(ax1,gcf);
        imagesc(ax1,Peak_Img/max(max(Peak_Img))); %% Plot the image
        imagesc(ax2,smoothdata2(tau_plot(:,:,1),"movmean",3), 'AlphaData', 0.7*(Mask > 0));
        clim(ax1,[param_plot.clim_I_min,param_plot.clim_I_max])
        % clim(ax2,[param_plot.clim_tau_min,param_plot.clim_tau_max])
        tau_temp2 = tau_plot(:,:,1);
        tau_temp3 = tau_temp2(Mask);
        tau_mean_temp = mean(tau_temp3(:));
        tau_std_temp = std(tau_temp3(:));
        clim(ax2,[max(0,tau_mean_temp-1.5*tau_std_temp),tau_mean_temp+1.5*tau_std_temp])
        colorbar(ax2)
        xlim([xmin ,xmax])
        ylim([ymin ,ymax])
        pbaspect(ax1,[xmax-xmin ymax-ymin 1])
        set(ax1, 'YDir', 'reverse');
        ax2.UserData = linkprop([ax1,ax2],{'Position','InnerPosition','DataAspectRatio','xtick','ytick','ydir','xdir','xlim','ylim'});
        ax3.UserData = linkprop([ax1,ax3],{'Position','InnerPosition','DataAspectRatio','xtick','ytick','ydir','xdir','xlim','ylim'});
        colormap(ax1,'gray')
        colormap(ax2,'turbo')
        
        if (param_plot.scale_Left_or_right == 1)
            x_vector = [xmin+20 xmin+Size_pixels+20];   % for plotting
        else
            x_vector = [xmax-20-Size_pixels xmax-20];   % for plotting
        end
        y_vector = [ymax-20 ymax-20];   % for plotting
        line(ax3,x_vector,y_vector,'LineWidth',param_plot.scale_bar_Linewidth,'Color',[0.99 0.99 0.99]);  % plot line
        set(ax1,'XColor', 'none','YColor','none')
        set(ax2,'XColor', 'none','YColor','none')
        set(ax3,'XColor', 'none','YColor','none')
        ax2.Visible = 'off';
        ax3.Visible = 'off';        
        fig_name = get_file_name('LTMap',fig_name_template,'.fig');
        png_name = get_file_name('LTMap',fig_name_template,param_plot.figformat);
    
        saveas(figure(4),fig_name)
        saveas(figure(4),png_name)
        close(figure(4));
    end
    
    %% Save colorbar separately
    if (param_plot.Save_colorbar == 1)
        gcf = figure(4);
        xmin = param_plot.Img_xlim_min;
        xmax = param_plot.Img_xlim_max;
        ymin = param_plot.Img_ylim_min;
        ymax = param_plot.Img_ylim_max;
    
        gcf.Position = [0,0,xmax-xmin,ymax-ymin];
        
        ax_temp = axes(gcf);hold on;
        ax_temp.InnerPosition = [0,0.1,0.8,0.8];
        colormap(ax_temp,'turbo')
        cb = colorbar(ax_temp);
        cb.Title.String = '\tau (ns)';
        cb.FontSize = 15 ;
        clim(ax_temp,[param_plot.clim_tau_min,param_plot.clim_tau_max])
        ax_temp .Visible = 'off';
        set(ax_temp ,'XColor', 'none','YColor','none')
        fig_name = strcat(fig_name_template,'LTMap_Cbar',param_plot.figformat);
        png_name = strcat(fig_name_template,'LTMap_Cbar.fig');
        saveas(figure(4),fig_name)
        saveas(figure(4),png_name)
        close(figure(4));
    end
    
    %% Save intensity images separately 
    if (param_plot.plot_PeakInt == 1)
        xmin = param_plot.Img_xlim_min;
        xmax = param_plot.Img_xlim_max;
        ymin = param_plot.Img_ylim_min;
        ymax = param_plot.Img_ylim_max;
        gcf = figure(param_plot.fig_idx);
        gcf.Position = [0,0,xmax-xmin,ymax-ymin];
        ax1 = axes(gcf);hold on;
        ax1.InnerPosition = [0,0,1,1];
        ax3 = copyobj(ax1,gcf);
        imagesc(ax1,Peak_Img/max(max(Peak_Img))); %% Plot the image
        clim(ax1,[param_plot.clim_I_min,param_plot.clim_I_max])
        xlim([xmin ,xmax])
        ylim([ymin ,ymax])
        pbaspect(ax1,[xmax-xmin ymax-ymin 1])
        
        set(ax1, 'YDir', 'reverse');
        ax3.UserData = linkprop([ax1,ax3],{'Position','InnerPosition','DataAspectRatio','xtick','ytick','ydir','xdir','xlim','ylim'});
        colormap(ax1,'gray')
        % colorbar;
        if (param_plot.scale_Left_or_right == 1)
            x_vector = [xmin+20 xmin+Size_pixels+20];   % for plotting
        else
            x_vector = [xmax-20-Size_pixels xmax-20];   % for plotting
        end
        y_vector = [ymax-20 ymax-20];   % for plotting
        line(ax3,x_vector,y_vector,'LineWidth',param_plot.scale_bar_Linewidth,'Color',[0.99 0.99 0.99]);  % plot line
        set(ax1,'XColor', 'none','YColor','none')
        set(ax3,'XColor', 'none','YColor','none')
        ax3.Visible = 'off';
        fig_name = get_file_name('IMap',fig_name_template,'.fig');
        png_name = get_file_name('IMap',fig_name_template,param_plot.figformat);
    
        saveas(figure(param_plot.fig_idx),fig_name)
        saveas(figure(param_plot.fig_idx),png_name)
        close(figure(param_plot.fig_idx));


        gcf = figure(4);
        xmin = param_plot.Img_xlim_min;
        xmax = param_plot.Img_xlim_max;
        ymin = param_plot.Img_ylim_min;
        ymax = param_plot.Img_ylim_max;
    
        gcf.Position = [0,0,xmax-xmin,ymax-ymin];
        
        ax_temp = axes(gcf);hold on;
        ax_temp.InnerPosition = [0,0.1,0.8,0.8];
        colormap(ax_temp,'gray')
        cb = colorbar(ax_temp);
        cb.Title.String = 'Intensity (counts)';
        cb.FontSize = 15 ;
        clim(ax_temp,[param_plot.clim_I_min*min(Peak_Img(:)),param_plot.clim_I_max*max(Peak_Img(:))])
        ax_temp .Visible = 'off';
        set(ax_temp ,'XColor', 'none','YColor','none')
        fig_name = strcat(fig_name_template,'LTMap_Cbar',param_plot.figformat);
        png_name = strcat(fig_name_template,'LTMap_Cbar.fig');
        saveas(figure(4),fig_name)
        saveas(figure(4),png_name)
        close(figure(4));
    end
    
    
    %% Automatically label and extract statistics
    if (param_plot.plot_AutoObjStats == 1)
        Data_to_write = cell(N_Object_auto,7);
        VariableNamesArray = cell(5+param.order*2,1);
        VariableNamesArray{1} = 'File Name';
        VariableNamesArray{2} = 'Object Index';
        VariableNamesArray{3} = 'Number of pixels';
    
        for idx_Obj = 1:N_Object_auto   %%  for each label separately
            str_temp = string(Get_savefile_name(param,path,"filename"));
            Data_to_write{idx_Obj,1} = str_temp;
            str_temp = char(str_temp);
            matches = regexp(str_temp, '_W(\d)', 'tokens');  % find the number after W
            if ~isempty(matches)
                Well_current = str2double(matches{1}{1});  % convert to number
            else
                Well_current = 0;
            end
            Mask_Temp = Obj_Mask;
            Mask_Temp(Mask_Temp ~= idx_Obj) = 0;
            Mask_Temp(Mask_Temp == idx_Obj) = 1;
            N_pixels_in_Obj = sum(Mask_Temp(:));
            Data_to_write{idx_Obj,2} = idx_Obj;
            Data_to_write{idx_Obj,3} = N_pixels_in_Obj;
            if ismember(param.Analysis_type,[1,2])
                %%%  Objects Histogram
                tau_Obj = tau_plot;
                if param.order == 1
                    tau_Obj(~Mask_Temp) = NaN;
                    [V,E] = histcounts(tau_Obj,param_plot.hbins_tau);
                    W = V/sum(V);
                    dE = diff(E(1:2));
                    tau_mean = sum((E(1:end-1)+dE/2).*W);
                    tau_std = sqrt(sum(W.*(E(1:end-1) + dE/2 - tau_mean).^2));
                    chi_Obj = chi_map;
                    chi_Obj(~Mask_Temp) = NaN;    
                    [V_chi,E_chi] = histcounts(chi_Obj ,param_plot.hbins_chi);
                    W_chi = V_chi/sum(V_chi); 
                    dE_chi = diff(E_chi(1:2));
                    chi_mean = sum((E_chi(1:end-1)+dE_chi/2).*W_chi);
                    chi_std = sqrt(sum(W_chi.*(E_chi(1:end-1) + dE_chi/2 - chi_mean).^2));
                    if exist('Photon_Budget','var')
                        PB_current = sum(Photon_Budget(logical(Mask_Temp)))/N_pixels_in_Obj;
                    else
                        PB_current = 0;
                    end
                    Data_to_write{idx_Obj,4} = tau_mean;
                    Data_to_write{idx_Obj,5} = tau_std;
                    Data_to_write{idx_Obj,6} = chi_mean;
                    Data_to_write{idx_Obj,7} = chi_std;
                    Data_to_write{idx_Obj,8} = Well_current;
                    Data_to_write{idx_Obj,9} = PB_current;
                    VariableNamesArray{4} = 'LT mean (ns)';
                    VariableNamesArray{5} = 'LT std (ns)';
                    VariableNamesArray{6} = 'chi mean';
                    VariableNamesArray{7} = 'chi std';
                    VariableNamesArray{8} = 'Well';
                    VariableNamesArray{9} = 'Photon Budget/pixel';
                    cellsize_col = 10;
                    celltype_col = 11;
                elseif param.order > 1
                    for idx_order = 1:param.order
                        tau_order_temp = tau_Obj(:,:,idx_order);
                        tau_order_temp(~Mask_Temp) = NaN;
                        [V,E] = histcounts(tau_order_temp,param_plot.hbins_tau);
                        W = V/sum(V);
                        dE = diff(E(1:2));
                        tau_mean_obj = sum((E(1:end-1)+dE/2).*W);
                        tau_std_obj = sqrt(sum(W.*(E(1:end-1) + dE/2 - tau_mean_obj).^2));
                        Data_to_write{idx_Obj,idx_order+3} = tau_mean_obj;
                        Data_to_write{idx_Obj,idx_order+3+param.order} = tau_std_obj;
                        VariableNamesArray{idx_order+3} = strcat('LT mean ',num2str(idx_order),' (ns)');
                        VariableNamesArray{idx_order+3+param.order} = strcat('LT std ',num2str(idx_order),' (ns)');
                    end
                    chi_Obj = chi_map;
                    chi_Obj(~Mask_Temp) = NaN;    
                    [V_chi,E_chi] = histcounts(chi_Obj ,param_plot.hbins_chi);
                    W_chi = V_chi/sum(V_chi); 
                    dE_chi = diff(E_chi(1:2));
                    chi_mean = sum((E_chi(1:end-1)+dE_chi/2).*W_chi);
                    chi_std = sqrt(sum(W_chi.*(E_chi(1:end-1) + dE_chi/2 - chi_mean).^2));
        
                    Data_to_write{idx_Obj,4+2*param.order} = chi_mean;
                    Data_to_write{idx_Obj,5+2*param.order} = chi_std;
                end
                VariableNamesArray{4+2*param.order} = 'chi mean';   
                VariableNamesArray{5+2*param.order} = 'chi std';
                VariableNamesArray{6+2*param.order} = 'Well';
                VariableNamesArray{7+2*param.order} = 'Photon Budget/pixels';
                cellsize_col  = 8+2*param.order;
                celltype_col  = 9+2*param.order;

            elseif ismember(param.Analysis_type,3)
                VariableNamesArray{4} = 'Well';
                VariableNamesArray{5} = 'tauUA';
                VariableNamesArray{6} = 'tauBA';
                VariableNamesArray{7} = 'wBM';
                VariableNamesArray{8} = 'wUA';
                VariableNamesArray{9} = 'wBA';
                VariableNamesArray{10} = 'rho_b';
                VariableNamesArray{11} = 'rho_c';
                VariableNamesArray{12} = 'wUM';

                Data_to_write{idx_Obj,4} = Well_current;
                Data_to_write{idx_Obj,5} = theta_avg{idx_Obj}(1);
                Data_to_write{idx_Obj,6} = theta_avg{idx_Obj}(2);
                Data_to_write{idx_Obj,7} = theta_avg{idx_Obj}(3);
                Data_to_write{idx_Obj,8} = theta_avg{idx_Obj}(4);
                Data_to_write{idx_Obj,9} = theta_avg{idx_Obj}(5);
                Data_to_write{idx_Obj,10} = theta_avg{idx_Obj}(6);
                Data_to_write{idx_Obj,11} = theta_avg{idx_Obj}(7);
                Data_to_write{idx_Obj,12} = 1 - sum(theta_avg{idx_Obj}(3:5));
                cellsize_col  = 13;
            end

            if (param.Stain)
                if (Stain_flag)
                    VariableNamesArray{cellsize_col} = 'CellSize';
                    VariableNamesArray{celltype_col} = 'CellType';
                    Data_to_write{idx_Obj,cellsize_col} = cellsize_arr(idx_Obj);
                    Data_to_write{idx_Obj,celltype_col} = celltype_array(idx_Obj);
                end
            end
        end


        if(isfile(Excel_path))
            newRowTable = cell2table(Data_to_write);
            writetable(newRowTable, Excel_path, 'WriteMode', 'append', 'Sheet', 1);
        else
            initialTable = cell2table(Data_to_write, 'VariableNames', VariableNamesArray);
            writetable(initialTable, Excel_path);
        end
    end
    drawnow;
    
    if (param_plot.plot_Colocalize == 1)
        if (param.Analysis_type == 3)
            I_D = Peak_Img_arr{1};
            I_A = Peak_Img_arr{3};

            % thresholds (set to 0; change if you want to exclude background)
            tD = 0;
            tA = 0;
            
            % masked vectors
            D = I_D(Mask);
            A = I_A(Mask);
            
            posD = D > tD;
            posA = A > tA;
            coloc = posD & posA;
            
            % Manders coefficients
            M1 = sum(D(coloc)) / sum(D(posD));   % donor fraction overlapping acceptor
            M2 = sum(A(coloc)) / sum(A(posA));   % acceptor fraction overlapping donor
            
            % (optional) Manders overlap coefficient (MOC)
            MOC = sum(D(coloc).*A(coloc)) / sqrt(sum(D(posD).^2) * sum(A(posA).^2));
            
            % Simple normalized composite for visualization (robust to outliers)
            pD = prctile(D(posD), 99);
            pA = prctile(A(posA), 99);
            
            I_Dn = zeros(size(I_D));
            I_An = zeros(size(I_A));
            I_Dn(Mask) = min(I_D(Mask)/pD, 1);
            I_An(Mask) = min(I_A(Mask)/pA, 1);
            
            RGB = zeros([size(I_D), 3]);
            RGB(:,:,1) = I_An;   % red  = acceptor
            RGB(:,:,2) = I_Dn;   % green= donor
            
            % Plot
            figure('Name','FRET colocalization (Manders)');
            
            subplot(1,3,1);
            imshow(RGB);
            title('Composite in Mask (R=A, G=D)');
            
            subplot(1,3,2);
            imagesc(Mask .* min(I_Dn, I_An));
            axis image off;
            colorbar;
            title('Colocalization map (min norm)');
            
            subplot(1,3,3);
            scatter(D(posD | posA), A(posD | posA), 6, 'k', 'filled'); hold on;
            scatter(D(coloc),       A(coloc),       6, 'r', 'filled');
            xlabel('Donor intensity');
            ylabel('Acceptor intensity');
            grid on;
            title(sprintf('Manders: M1=%.3f, M2=%.3f | MOC=%.3f', M1, M2, MOC));
            
            % "Manders correlation" plot (pixelwise normalized product inside mask)
            % This is a visual correlation-like heatmap (not a separate coefficient).
            figure('Name','Manders correlation map');
            imagesc(Mask .* (I_Dn .* I_An));
            axis image off;
            colorbar;
            title(sprintf('I_Dn .* I_An inside Mask | M1=%.3f, M2=%.3f', M1, M2));


        end
        drawnow;
    end
    %% Filename function
    function fig_name = get_file_name(name_str,fig_name_template,str_formate)
        fig_name_temp = strcat(fig_name_template,name_str,str_formate);
        i = 1;
        while(isfile(fig_name_temp))
            fig_name_temp = strcat(fig_name_template,name_str,'_',num2str(i),str_formate);
            i = i+1;
        end
        fig_name = fig_name_temp;
    end
    
    function [histw, histv] = histwv(v, w, min, max, bins)
    
        %Inputs: 
        %vv - values
        %ww - weights
        %minV - minimum value
        %maxV - max value
        %bins - number of bins (inclusive)
        
        %Outputs:
        %histw - wieghted histogram
        %histv (optional) - histogram of values    
       
        delta = (max-min)/(bins-1);
        subs = round((v-min)/delta)+1;
        
        histw = accumarray(subs(:),w(:),[bins,1]);
        if nargout == 2
            histv = accumarray(subs(:),1,[bins,1]);
        end
            
    end
    
    function [chi_mean,chi_std,chi_dist] = FitHistogram_chi(x,W,chi_mean,chi_std)
        
        get_dist = @(b) normpdf(x,b(1),b(2))/sum(normpdf(x,b(1),b(2)));
        Costfn = @(b) sum(((get_dist(b) - W).^2)./1) + 1000000*(b(1)<0);
        IC = [chi_mean,chi_std];
        
        B_final  = fminsearch(Costfn,IC);
        chi_mean = B_final(1);
        chi_std = B_final(2);
        chi_dist  = get_dist(B_final);
    end
    
    function [I_amp,I_mu,I_std,exp_dist,x] = FitHistogram_tau(x, I, param,IC)
    
        W = (I/sum(I));
        W = reshape(W,[1,length(W)]);
        
        if (param.order == 1)
            get_dist = @(b) normpdf(x,b(1),b(2))/sum(normpdf(x,b(1),b(2)));
            Costfn = @(b) sum(((get_dist(b) - W).^2)./1) + 1000000*(b(1)<0);
    
        elseif (param.order == 2)
            get_dist = @(b) b(1)*normpdf(x,b(2),b(3))/sum(normpdf(x,b(2),b(3))) + (1-b(1))*normpdf(x,b(4),b(5))/sum(normpdf(x,b(4),b(5)));
            Costfn = @(b) sum(((get_dist(b) - W).^2)./1) + 1000000*(b(1)<0 || b(2) < 0 || b(4) < 0) + 1000000*(b(3)>0.8 || b(5) > 0.8 );
            Amp = IC(1)/(IC(1)+IC(4));
            IC = [Amp IC(2) IC(3) IC(5) IC(6)];
        end
        
        B_final  = fminsearch(Costfn,IC);
    
        if (param.order == 2 && (B_final(1) < 0.15 || B_final(1) > 0.85 || B_final(3) > 1 || B_final(5) > 1 || (B_final(1) < 0.2 && B_final(3) > 0.4) || (B_final(1) > 0.8 && B_final(5) > 0.4)  ))  % Check if the second component is negligible
            % If the weight of the second component is too small, switch to order 1
    
            get_dist = @(b) normpdf(x,b(1),b(2))/sum(normpdf(x,b(1),b(2)));
            Costfn = @(b) sum(((get_dist(b) - W).^2)./1) + 1000000*(b(1)<0);
            B_final  = fminsearch(Costfn,IC);
            I_amp = [0.5 0.5];
            I_mu = [B_final(1) B_final(1)];
            I_std = [B_final(2) B_final(2)];
            exp_dist  = get_dist(B_final);
    
            % if (B_final(1) > 0.85 || B_final(5)>1)
            %     I_mu  = [B_final(2) B_final(2)];
            %     I_std = [B_final(3) B_final(3)];
            % elseif (B_final(1) < 0.15 || B_final(3)>1)
            %     I_mu  = [B_final(4) B_final(4)];
            %     I_std = [B_final(5) B_final(5)];
            % end
    
        else
            if (param.order == 1)
                I_amp = 1;
                I_mu = B_final(1);
                I_std = B_final(2);
            elseif (param.order == 2)
                I_mu  = [B_final(2) B_final(4)];
                [I_mu idx_sort] = sort(I_mu);
                I_amp = [B_final(1) (1-B_final(1))];
                I_std = [B_final(3) B_final(5)];
                I_amp = I_amp(idx_sort);
                I_std = I_std(idx_sort);
            end
            exp_dist  = get_dist(B_final);
        end
     end
    
    function CRB_out = GetCRB(tau,param)
    
        amplist = [0.5,0.5];
        %%%% tau = [tau1 tau2]; %%%  We get the CRB of tau1
        amp = 2500;
        path_IRF = '/home/min/a/salem8/lib/Fitting/IRF50t_200G_2BIN_SP_1D.mat';
        load(path_IRF,"IRF","t_IRF");
        common_time = logical([repmat([1,0,0,0],1,125) 1])';
        t_IRF = t_IRF';
        t_sig = t_IRF(common_time);
        H = IRF;
        H = H-mean(H(t_IRF-t_IRF(1)<3.2));
        H = H/sum(H(:));
        t = t_IRF;
        t = t - t(1);
    
        B = [0 amp];
        for k = 1:param.order
            B(3+(k-1)*2) = tau(k);
            if (k < param.order)
                B(4+(k-1)*2) = amplist(k);
            end
        end
        B = [B 210];
        %%%% Salem: Take care of this part. It's fixed at 600
        MCP = 600;
        load(strcat('Noise_histogram_take5_',num2str(MCP),'.mat'),'Int_var_arr')
        y_current = output_t(H,t,B,param,common_time);
    
        %%% Construct partial derivatives
        for i = 1:length(B)
            dy(:,i) =  y_deriv(H,t,B,param,i,common_time);
        end
    
        %%% Constructing the covariance matrix, should be diagonal if noise is independent
        Cm = zeros(size(y_current,1));
        for j = 1:size(y_current,1)
            Cm(j,j) = Int_var_arr(cast(abs(y_current(j)),"uint16"))';
        end
    
        %%% Computing Fisher information matrix
        I = zeros(length(B));    
        for i = 1:length(B)
            for j= 1:length(B)
                I(i,j) = transpose(dy(:,i))*inv(Cm)*dy(:,j);
            end
        end
        CRB = sqrt(inv(I));
        CRB_out = CRB(3,3);
    
        %% Function for forward model
        function dydp =  y_deriv(H,t,B,param,param_idx,common_time)
        %%% Should the dp step depend on the current value of the parameter?
            y1 = output_t(H,t,B,param,common_time);
            dp = B(param_idx)+0.01;
            B(param_idx)  = B(param_idx)+0.01;
            y2 = output_t(H,t,B,param,common_time);
            dydp = (y2-y1)/0.01;
        end
        
        function output = output_t(h,t,B,param,common_time)
            %%% Parameters are: [shift Amp  mu   sigma amp  ... ]
            %%% Parameters are: [b(1)  b(2) b(3) b(4)  b(5) ...]
            dt = diff(t(1:2));
            exp_vector = zeros(size(h,1),1);
        
            if (param.cont == 1)
                exp_dist = get_dist(B);
                exp_vector = sum(exp_dist .* exp(-t./ param.tau_vec),2);
            elseif (param.cont == 0)
        
                amp_temp = 0;
                for k = 1:param.order
                    if (k == param.order)
                        amp = 1-amp_temp;
                    else
                        amp = B(4+(k-1)*2);
                        amp_temp = amp_temp + amp;
                    end
                    exp_vector = exp_vector + amp*exp(-t/B(3+(k-1)*2));
                end
            end
        
        
            exp_vector(isnan(exp_vector)) = 0;    %%% To avoid any errors when first point is NaN at 0 lifetime
            output  = myconv(h,exp_vector,B(1));
            output  = B(2)*output/max(output)+B(end);
            output  = output(common_time);
        
            function exp_dist = get_dist(B)
                tau_vector = param.tau_vec;
                exp_dist = zeros(size(tau_vector));
                if (param.cont == 0)
                    for j = 1:param.order
                        tau_temp = B(3+2*(j-1));
                        idx = find(abs(tau_vector - tau_temp) == min(abs(tau_vector - tau_temp)));
                        exp_dist(idx) = exp_dist(idx) + B(4+(j-1)*2);
                    end
                else
                    for j = 1:param.order
                        mu = B(3+(j-1)*3);
                        sigma = B(4+(j-1)*3);
                        if (sigma == 0)
                            continue;
                        end
                        if(param.cont_type == 1)    %%%% Gaussian distribution
                            pdf_temp = normpdf(tau_vector,mu,sigma)*(tau_vector(end)-tau_vector(1))/length(tau_vector);
                        elseif(param.cont_type == 2)    %%%% Gama distribution
                            pdf_temp = gampdf(tau_vector,mu^2/sigma,sigma/mu)*(tau_vector(end)-tau_vector(1))/length(tau_vector);
                        elseif(param.cont_type == 3)    %%%% Inverse Gama distribution
                            pdf_temp = gampdf(1./tau_vector, mu^2 / sigma, sigma/(mu^3)) * (tau_vector(end) - tau_vector(1)) / length(tau_vector);    
                        elseif(param.cont_type == 4)    %%%% Log-normal distribution
                            gamma = asinh(sigma/(2*mu));
                            pdf_temp = exp(-((log(1./tau_vector)-log(1/mu))/gamma).^2);    
                        end
                        exp_dist = exp_dist + B(5+(j-1)*3)*pdf_temp;
                    end
                end
            end
            
            function C = myconv(A, B,t_shift)   % perform 1D Linear discrete convolution followed by clipping
                N = length(A);
                %%% Zero Padding
                N_ZP = 1024;
                A_ZP = cat(1,A,zeros(N_ZP-N,1));
                B_ZP = cat(1,B,zeros(N_ZP-N,1));
                f_max =1/dt;
                df_ZP=1/(N_ZP*dt);
                f_ZP = (-f_max/2:df_ZP:f_max/2-df_ZP)';
                %%% FFT --> shift and Conv --> IFFT
                A_W = FT(A_ZP,dt);
                B_W = FT(B_ZP,dt);
                try
                    C_W = A_W.*B_W.*exp(-1j*2*pi*f_ZP*t_shift);
                catch
                    if (param.DEBUGGING == 1)
                        keyboard;
                    else
                        err_status = 1;
                        return;
                    end
                end
                C = abs(IFT(C_W,1/(N*dt)));
                %%% Clipping the zero padded signal
                C = C(1:N);
            end
            
            function G = FT(g, dt)  % 1D DFT where delta is sample width
                G = fftshift(fft(fftshift(g))) * dt;
            end
            
            function g = IFT(G, df) % 1D DIFT where df is sample with in frequency domain (1/N/dt)
                g = (ifft(ifftshift(G)) * length(G) * df);    
            end
            
        end 
    end


end

%% Other functions

function str = Get_savefile_name(param,path,mode)
    filename_saving_temp = split(strip((path)),'/');
    current_foldername = strcat(filename_saving_temp(end-2),"_",filename_saving_temp(end-1));
    current_filename = filename_saving_temp(end);
    folder1 = fullfile('Results',current_foldername);
    folder2 = fullfile('Results',current_foldername,current_filename);
    if (~isfolder(folder1))
        mkdir(char(folder1));
    end
    if (~isfolder(folder2))
        mkdir(char(folder2));
    end

    if (param.cont == 1)
        str_cont = strcat('C',num2str(param.cont_type));
    else
        str_cont = 'D';
    end
    if (param.SolveMatrix == 1)
        MatrixStr = "M";
    else
        MatrixStr = [];
    end

    if (param.fixed_amp == 1 || param.fixed_amp_prior == 1 && ~strcmp(mode,"Results1D"))
        str_Amp_prior = 'A';
    else
        str_Amp_prior = [];
    end

    if (param.Use_Gradient == 1)
        str_Gradient = 'G';
    else
        str_Gradient = [];
    end

    if (param.EM_Algorithm == 1)
        str_EM = 'EM';
    else
        str_EM = [];
    end
    if (~isempty(param.NameSTR))
        str_N = strcat('_',param.NameSTR);
    else
        str_N = [];
    end
    %%% Setting mask settings
    if (param.Analysis_type == 2)
        if (param.PriorMask(param.folder_idx_current) > 0)
            str_PriorMask = "U";
        elseif (param.PriorMask(param.folder_idx_current) < 0)
            str_PriorMask = "S";
        else
            str_PriorMask = [];
        end
    else
        str_PriorMask = [];
    end
    %%% Setting Prior settings
    if (~isempty(param.PriorCell) || ~isempty(param.PriorData) &&  ~strcmp(mode,"Results1D"))
        str_prior = 'P';
    elseif (param.Analysis_type == 2 &&  ~strcmp(mode,"Results1D"))
        if (param.PriorVector(param.folder_idx_current) > 0)
            str_prior = 'P';
        else
            str_prior = [];
        end
    else
        str_prior = [];
    end


    if strcmp(mode,"Results")
        if ismember(param.Analysis_type,[1,2])
            str = strcat('Results/RESULTS_EXP',num2str(param.order),str_cont,'_',current_foldername,current_filename,'_M',num2str(param.Method),str_Gradient,str_prior,str_PriorMask,str_Amp_prior,MatrixStr,str_N,str_EM,'.mat');
        elseif ismember(param.Analysis_type,3)
            str = strcat('Results/RESULTS_Cell_EXP',num2str(param.order),str_cont,'_',current_foldername,current_filename,'_M',num2str(param.Method),str_Gradient,str_prior,str_PriorMask,str_Amp_prior,MatrixStr,str_N,str_EM,'.mat');
        end
    elseif strcmp(mode,"Results1D")
        str = strcat('Results/RESULTS_EXP1',str_cont,'_',current_foldername,current_filename,'_M',num2str(param.Method),str_Gradient,str_prior,str_PriorMask,str_Amp_prior,MatrixStr,str_N,str_EM,'.mat');
    elseif strcmp(mode,"Mask")
        str = strcat('Results/Masks/MaskThresh_',current_foldername,'_',current_filename,str_PriorMask,'.mat');
    elseif strcmp(mode,"AmpMap")
        str = strcat('Results/AmpMap/AmpMap_',current_foldername,'_',current_filename,str_PriorMask,'.mat');
    elseif strcmp(mode,"Label")
        str = strcat('Results/Masks/MaskLabel_',current_foldername,'_',current_filename,str_PriorMask,'.mat');
    elseif strcmp(mode,"Excel")
        if ismember(param.Analysis_type,[1,2])
            str = strcat('Results/',current_foldername,'/Stats_',current_foldername,'_Order',num2str(param.order),str_prior,str_EM,'.xlsx');
        elseif ismember(param.Analysis_type,3)
            str = strcat('Results/',current_foldername,'/Stats_Cell_',current_foldername,'_Order',num2str(param.order),str_prior,str_EM,'.xlsx');
        end
    elseif strcmp(mode,"Template")
        str = strcat(folder2,'/','M',num2str(param.Method),str_Gradient,str_prior,str_PriorMask,str_Amp_prior,MatrixStr,'_exp',num2str(param.order),str_cont);
    elseif strcmp(mode,"filename")
        str = current_filename;
    else
        str = [];
    end
    str = string(str);
end


function [C,str_leg,str_linestyle,str_MarkerStyle] = Get_Method_info(param)

    Method = param.Method;
    if (~isempty(param.PriorCell) || ~isempty(param.PriorData))
        if (param.fixed_amp == 1 || param.fixed_amp_prior == 1)
            Method = 18;
        else
            Method = 17;
        end
    end
    str_linestyle = '-';
    str_MarkerStyle = 'x';
    switch Method
        case 2
            C = [1,0,1];
            str_leg = "Gaussian with shot noise";
        case 3
            C = [0,85,212]/255;
            str_leg = "Poisson Model";
        case 6
            C = [0.47,0.67,0.19];
            str_leg = "Empirical Model";
        case 7
            C = [0.72,0.27,1.00];
            str_leg = "Neighbors Model";
        case 8
            C = [0,0,0];
            str_leg = "True Noise";
            str_linestyle  = '--';
        case 9
            if (~isempty(param.PriorCell))
                C = [0,1,0.5];
                str_leg = "Hybrid Model + Prior";
            else
                C = [255,105,41]/255;
                str_leg = "Hybrid Model";
                str_MarkerStyle = 'o';
            end
        case 10
            C = [0,0,0];
            str_leg = "True Noise";
        case 11
            if (~isempty(param.PriorCell))
                C = [120,171,48]/255;
                str_leg = "Hybrid Model + Taper + Prior";
                str_MarkerStyle = 'square';
            else
                C = [125,46,143]/255;
                str_leg = "Hybrid Model + Taper ";
                str_MarkerStyle = '*';
            end
        case 13
            C = [0,1,0.5];
            str_leg = "Hybrid Model + Prior";
        case 17
            C = [120,171,48]/255;
            str_leg = "Hybrid Model + Prior";
        case 18
            C = [125,46,143]/255;
            % C = [0.07,0.62,1.00];
            str_leg = "Hybrid Model + Prior";

        case 20
            C = [0,0,1];
            % C = [0.07,0.62,1.00];
            str_leg = "Hybrid Model + Prior";
    end
    % if (param.order == 1)
    %     C = [0,0,0];
    % end
end



function [y_fit,Noise_var,amp,tau,chi_out,err_status,chi_vec,Ar_out,DC_out,tshift_out,err_str,param] ... 
    = mydeconv8(y ,w, t_sig, H, t_IRF,param,NoiseVector,y_neighbors,w_neighbors,y_avg_tshift,Prior_Current)
    
    if (~isempty(param.PriorCell) && param.randPrior)
        rand_arr = rand(2*param.order,1);
    end
    t_sig = double(cast(t_sig*1000,'int64'));
    t_IRF = double(cast(t_IRF*1000,'int64'));
    %%% Defining zero vectors for outputs
    Noise_var = [];
    tau = zeros(param.order,1);
    amp = zeros(param.order,1);
    sig = zeros(param.order,1);
    err_status = 0;
    err_str = [];
    chi_out = [];
    chi_vec = [];
    Ar_out = [];
    DC_out = [];
    y_fit = [];
    tshift_out = [];
    if (param.Method == 20)
        curve_fitting_flag = 0;
        phasor_flag = 1;
    else
        curve_fitting_flag = 1;
        phasor_flag = 0;
    end

%%% Defining some parameters
    fontsz = 9;
    dt_sig = diff(t_sig(1:2));
    dt_IRF = diff(t_IRF(1:2));
    N = length(y); 
    f_max = 1000/dt_IRF; %%% the 1000 will account for the change from ns to ps
    df= 1000.0/(N*dt_IRF*1.0);
    f = (-f_max/2:df:f_max/2-df)';
    N_ZP = 1024;
    df_ZP = 1000/(N_ZP*dt_IRF);
    f_ZP = (-f_max/2:df_ZP:f_max/2-df_ZP)';
    % f_ZP = (-N_ZP/2 : N_ZP/2 - 1)' * df_ZP;  % ← use this exact form

%%% Rescaling
    if (param.laser == "SP")
        H = H-mean(H(t_IRF-t_IRF(1)<1000*3.2));
        % DarkCurrent = mean(y(t_sig-t_IRF(1)<1000*3.3));
    elseif (param.laser == "NKT")
        H = H-mean(H(t_IRF-t_IRF(1)<1000*3.2));
        % DarkCurrent = mean(y(t_sig-t_IRF(1)<1000*3.2));
    elseif (param.laser == "OldPeak30")
        H = H-mean(H(t_IRF-t_IRF(1)<1000*2.5));
    elseif (param.laser == "OldPeak6")
        H = H-mean(H(t_IRF-t_IRF(1)<1000*0.6));
    elseif (param.laser == "NKTOldPeak22")
       % H = H-mean(H(t_IRF-t_IRF(1)<1000*1));
    end
    DarkCurrent = param.DCShift;
    if (param.SubDrk == 1 && param.Method ~= 2)
        y = y + DarkCurrent;
        y_neighbors = y_neighbors + DarkCurrent;
    elseif (param.SubDrk == 1 && param.Method == 2) 
        DarkCurrent = 5;
        y = y + DarkCurrent;
        y_neighbors = y_neighbors + DarkCurrent;        
    end
    H = H/sum(H);
    if (param.Method == 3)
        y(y<1) = 1;
    end

    idx_peak = find(abs(y - max(y)) == 0);
    if (param.after_peak == 0)
        idx_peak = -2;
    end

    if (param.clip_type == 1 && param.Method == 9)
        diff_arr = FindDiff(y,1);
        ratio_arr = sqrt(w).*(y-diff_arr);
        clip_vector = (ratio_arr < 4);
    elseif (param.clip_type == -1 && param.Method == 9)
        y_temp = smooth(y,5);
        clip_vector = (y_temp > 350);
    else
        clip_vector = idx_peak+3:1:length(y);
    end

%% After this point, both the signal in time axes are all aligned
    t_sig = t_sig-t_sig(1); %%% Used in minimization
    t_IRF = t_IRF-t_IRF(1); %%% Used in forward model
    %%% Multiplication by 1000 helps stabilize the next line by avoiding
    %%% float points comparison
    common_time = ismember(t_IRF,t_sig);
    t_sig = t_sig/1000;
    t_IRF = t_IRF/1000;
    t = t_sig;
    dt_IRF = dt_IRF/1000;
    dt_sig = dt_sig/1000;
    % NoiseVector(NoiseVector<5) = 5;
%%% Calculating Scaling parameters using a mono exponential
    Ar_const = 1;
    if param.NormDivTau == 1
        temp_e = myconv(H,exp(-t_IRF/param.IC)/param.IC,0);     %%% Using tau = 4 ns
    else
        temp_e = myconv(H,exp(-t_IRF/param.IC),0);     %%% Using tau = 4 ns
    end
    [Ar_IC,tshift_temp] = align_rising_edge(y_avg_tshift,temp_e,DarkCurrent);
    tshift = tshift_temp;
%%
    if (param.Method == 8 || param.Method == 20)
        B_Actual = [];
        if (param.tshift_flag == 1)
            B_Actual = [B_Actual 0];
        end
        if (param.Scaleshift_flag == 1)
            B_Actual = [B_Actual Ar_IC];
        end
        for idx_order_k = 1:param.order
            if (param.cont == 0)
                if (param.order == 1)
                    B_Actual = [B_Actual param.NoiseActual_params(idx_order_k)];
                else
                    B_Actual = [B_Actual param.NoiseActual_params(idx_order_k) 1/param.order];
                end
            else
                if (param.order == 1)
                    B_Actual = [B_Actual param.NoiseActual_params(idx_order_k) 0.01];
                else
                    B_Actual = [B_Actual param.NoiseActual_params(idx_order_k) 0.01 1/param.order];
                end
            end
        end
        if (param.DCshift_flag == 1)
            B_Actual = [B_Actual 0];
        end
        % tshift = 0;
        y_Actual_temp = output_t(B_Actual);
        y_Actual = (param.Amp_Actual-DarkCurrent)*(y_Actual_temp-DarkCurrent)/max((y_Actual_temp-DarkCurrent))+DarkCurrent;
    end
    tshift = tshift_temp;
    y_OG = y;    
    idx_peak_chiout = idx_peak;
    y_all = [y_neighbors;y_OG']';
    N_all = size(y_all,2);
    NoiseVector_empirical  = param.MCPNoise(cast(abs(y_OG),"int16"));
    if (param.Averaging == 0)
        NoiseVector_empirical_neighbors = 1./w; 
    else
        NoiseVector_empirical_neighbors = sum(1./w_neighbors)/(N_all); 
        y_OG = mean(y_all,2);    
        y = y_OG;
    end
    if (param.Method ~= 8 && param.DEBUGGING == 0)
        y_Actual = y_OG; 
    end
    if (param.Method == 8)
        NoiseActual = param.MCPNoise(cast(abs(y_Actual),"int16"));
    else
        NoiseActual = NoiseVector_empirical_neighbors;
    end


    %% For debugging
    % Debug_plots();
    % if (param.DEBUGGING == 1 && param.PLOT_FLAG == 1)
    %     figure(111);hold on;
    %     plot(t_sig,y_OG./max(y_OG),'r','linewidth',2)
    %     plot(t_IRF,temp_e/max(temp_e),'b','linewidth',2)
    % end

    %% Clip noise 
    
    NoiseVector_empirical_clipped   = NoiseVector_empirical(clip_vector);
    NoiseVector_empirical_neighbors_clipped = NoiseVector_empirical_neighbors(clip_vector);
    if (param.NeighborsNoise_flag == 1 || param.Method == 7 || param.Method == 8)
        NoiseVector_clipped  = NoiseVector(clip_vector);
        NoiseActual_clipped = NoiseActual(clip_vector);
    end
    t = t(clip_vector);
    y_OG = y_OG(clip_vector);
    y = y_OG;
%% Defining the problem
if (curve_fitting_flag) %%% Curve fitting
    [problem,IC] = get_problem();   %%% IC here for debugging
    if (param.global == 1)
        rng default % For reproducibility
        gs = GlobalSearch;
        [B,Cost]   = run(gs,problem);
    else
        try
            [B,Cost] = fmincon(problem);
        catch ME
            if strcmp(ME.identifier, 'MATLAB:nearlySingularMatrix')
                disp('⚠️ Warning: Nearly singular matrix detected!');
                keyboard  % pause here
            else
                rethrow(ME)
            end
        end
    end
    y_fit  = output_t(B);
    Ar_out = Get_param(B,param,"Ar");
    DC_out = DarkCurrent;
    tshift_out = Get_param(B,param,"tshift");

%% Calculating chi-square of the fitted result
    if (param.Method == 2)
        chi_vec = ((y_fit - y_OG).^2)./param.MCPNoise(cast(abs(y_fit),"int16"));
    elseif (param.Method == 3 || param.Method == 12)  %%% Poission distribution
        chi_vec  = 2*(y_OG.*log(y_OG./y_fit)-(y_OG-y_fit));
        Noise_var = y_OG;
    elseif (param.Method == 5)
        chi_vec = ((y_fit - y_OG).^2)./(B(end-1)*y_OG);
    elseif (param.Method == 6)
        chi_vec = ((y_fit - y_OG).^2)./NoiseVector_empirical_clipped;
        Noise_var = NoiseVector_empirical_clipped;
    elseif (param.Method == 7)
        chi_vec = ((y_fit - y_OG).^2)./NoiseVector_clipped;
    elseif (param.Method == 8 || param.Method == 10)
        chi_vec = ((y_fit - y_OG).^2)./NoiseActual_clipped;
    elseif (param.Method == 9 || param.Method == 11 || param.Method == 13)
        chi_vec = ((y_fit - y_OG).^2)./NoiseVector_empirical_neighbors_clipped;
        Noise_var = NoiseVector_empirical_neighbors_clipped;
    elseif (param.Method == 14)
        m = param.NoiseCoeff;
        v = param.NoiseCoeff;
        mu_pearson = param.NoiseCoeff;
        sigma_pearson = param.NoiseCoeff;
        chi_vec  = log(abs((gamma(m+1i*v/2))./(gamma(m)))./(sigma_pearson*beta(m-0.5,0.5)))-m*log(1-((y-mu_pearson)/sigma_pearson).^2)-v*atan((y-mu_pearson)/sigma_pearson);
    end
    
    chi_r_vec = chi_vec/(length(chi_vec)-length(B));
    chi = sum(abs(chi_vec));
    chi_r = sum(abs(chi_r_vec));
    chi_vec_clipped = chi_vec(real(idx_peak_chiout)+3:end);
    chi_r_clipped = sum(chi_vec_clipped)/(length(chi_vec_clipped)-length(B));
    chi_out = chi_r_clipped;
    param.tau_dist = get_dist(B);

%% Fill output vectors with fitting results
    for idx_k = 1:param.order
        tau(idx_k) = Get_param(B,param,"tau",idx_k);
        amp(idx_k,1) = Get_param(B,param,"amp",idx_k);
        if (param.cont == 1)
            sig(idx_k) = Get_param(B,param,"sig",idx_k);
        end
    end

    if (param.SolveMatrix == 1)
        [~,A] = output_t(B);
        for k1 = 1:param.order
            tau(k1) = Get_param(B,param,"tau",k1);
            amp(k1) = A(k1+1);
        end
    end
    %% Print results and plot if required flags are on
    if (param.DEBUGGING == 1)
        disp(B);
        disp(strcat('Cost: ' , num2str(Cost) , '. chi:', num2str(chi),'. chi_r: ', num2str(chi_r),'. chi_peak',num2str(chi_out)));
        disp(strcat('Lifetime: ',num2str(tau),', amplitude: ',num2str(amp),', Sigma: ',num2str(sig)));
        if(param.PLOT_FLAG == 1)
            plot_I;
            keyboard    %%% Pauses to see the result before moving to the next one
        end
    end

elseif (phasor_flag)
    tau = phasor_analysis(y, t,1e6,param.NoiseActual_params(1));
end
    % [tau, idx ] = sort(tau);  %%% Sort B for better viewing
    % amp = amp(idx);
    
%% Functions used
    function varargout = Costfn(b)
        if param.Use_Gradient == 1
            grad_arr = zeros(length(b),1);  % Match your sampling axis size
        end
        if (param.SolveMatrix == 1)
            [y_calc,A] = output_t(b);
        else
            if param.Use_Gradient == 1
                [y_calc,grad_y_arr] = output_t(b);
            else
                y_calc = output_t(b);
            end
        end
        % y_calc = y_calc(clip_vector);
        y_meas = y;%(clip_vector);
        
        if (param.N_points ~= 0)
            %%% In here we start from the peak and measure N points
            %%% afterwards
            t_clip = t(clip_vector);
            if (isempty(param.N_points_arr))
                param.N_points_arr = 1:param.N_points;                
            end
            sampling_arr = reshape((param.N_points_arr+t_clip(1)),[1,length(param.N_points_arr)]);  %% Reshape used to gurantee dimensions are correct
            diff_matrix = abs(t_clip-sampling_arr);
            [~, idx_sampled] = min(diff_matrix, [], 1);
            y_calc = y_calc(idx_sampled);
            y_meas = y_meas(idx_sampled);
        else
            idx_sampled = 1:length(y_meas);
        end
        % 
        % if (param.Use_Gradient == 1 && param.tshift_flag == 1)
        %     delta = 1e-6;
        %     B_plus = b;
        %     B_minus = b;
        %     B_plus(1) = b(1) + delta;
        %     B_minus(1) = b(1) - delta;
        %     y_plus = output_t(B_plus);
        %     y_minus = output_t(B_minus);
        %     grad_y_arr(:,1) = (y_plus - y_minus) / (2*delta);
        % 
        % end
        % 
        % if param.DEBUGGING == 1 && param.Use_Gradient == 1
        %     delta = 1e-6;
        %     grad_numeric = zeros(size(b));
        %     for i = 1:length(b)
        %         B_plus = b;
        %         B_minus = b;
        %         B_plus(i) = b(i) + delta;
        %         B_minus(i) = b(i) - delta;
        % 
        %         y_plus = output_t(B_plus);
        %         y_minus = output_t(B_minus);
        % 
        %         grad_numeric(i) = sum(y_plus - y_minus) / (2*delta);
        %             figure(1);clf;hold on;
        %             plot(t_sig,grad_y_arr(:,i),'Linewidth',2)
        %             plot(t_sig,(y_plus - y_minus) / (2*delta),'r--','LineWidth',2)
        %             legend('Analytical','Numerical'); title('tshift Gradient Check');
        %             keyboard
        %     end
        %     grad_y_arr_temp = sum(grad_y_arr,1);
        %     fprintf("Parameter Gradient Comparison:\n");
        %     disp(table((1:length(b))', grad_numeric(:), grad_y_arr_temp(:), ...
        %         abs(grad_numeric(:) - grad_y_arr_temp(:)), ...
        %         'VariableNames', {'ParamIdx','Numerical','Analytical','AbsError'}));
        % 
        %     fprintf("Max abs diff: %g\n", max(abs(grad_numeric(:) - grad_y_arr_temp(:))));
        %     keyboard;
        % end



        switch param.Method
        case 1  %%% LMS
            C = sum((y_calc - y_meas).^2);
            if (param.Use_Gradient == 1)
                grad_arr = 2*sum(grad_y_arr.*(y_calc-y_meas));
            end
        case 2  %%% Gaussian noise, weighted LMS, dark current separated
            if max(y_calc) > 4095
                var_arr = interp1(1:1:4095,param.MCPNoise,1:1:ceil(max(y_calc)),"linear","extrap");
                var_y = var_arr(cast(abs(y_calc),"int16"))';
            else
                var_y = param.MCPNoise(cast(abs(y_calc),"int16"));
            end
            C = sum(((y_calc - y_meas).^2)./var_y + log(var_y)) + log(2*pi);
            if (param.Use_Gradient == 1)
                delta = 1;
                var_arr = interp1(1:1:4095,param.MCPNoise,1:1:ceil(max(y_calc))+1,"linear","extrap")';
                f_plus = var_arr(cast(abs(y_calc + delta), "int16"));
                f_minus = var_arr(cast(abs(y_calc - delta), "int16"));
                dvar = (f_plus - f_minus) / (2 * delta);
                
                dC_dy = (2 * (y_calc - y_meas) ./ var_y) ...
                      - (((y_calc - y_meas).^2) ./ var_y.^2) .* dvar ...
                      + (1 ./ var_y) .* dvar;

                grad_arr = sum(grad_y_arr.*dC_dy);
            end
        case {3,12}  %%% MLE with poisson prior distribution
            C = 2*sum(y_meas.*log(y_meas./y_calc)-(y_meas-y_calc));
            if (param.Use_Gradient == 1)
                grad_arr = 2*sum(grad_y_arr.*(1-y_meas./y_calc));
            end
        case 6  %%% Gaussian noise, weighted LMS
            C = sum(abs(((y_calc - y_meas).^2)./NoiseVector_empirical_clipped(idx_sampled)));
            if (param.Use_Gradient == 1)
                grad_arr = 2*sum(grad_y_arr.*(y_calc-y_meas)./NoiseVector_empirical_clipped(idx_sampled));
            end
        case 7
            C = sum(abs(((y_calc - y_meas).^2)./NoiseVector_clipped(idx_sampled)));
            if (param.Use_Gradient == 1)
                grad_arr = 2*sum(grad_y_arr.*(y_calc-y_meas)./NoiseVector_clipped(idx_sampled));
            end

        case {8,10}
            C = sum(abs(((y_calc - y_meas).^2)./NoiseActual_clipped(idx_sampled)));
            if (param.Use_Gradient == 1)
                grad_arr = 2*sum(grad_y_arr.*(y_calc-y_meas)./NoiseActual_clipped(idx_sampled));
            end
        case {9,11,13}
            C = sum(abs(((y_calc - y_meas).^2)./NoiseVector_empirical_neighbors_clipped(idx_sampled)));
            if (param.Use_Gradient == 1)
                grad_arr = 2*sum(grad_y_arr.*(y_calc-y_meas)./NoiseVector_empirical_neighbors_clipped(idx_sampled));
            end
        case 14
            m = param.NoiseCoeff;
            v = param.NoiseCoeff;
            mu_pearson = param.NoiseCoeff;
            sigma_pearson = param.NoiseCoeff;
            C = sum(-m*log(1-((y-mu_pearson)/sigma_pearson).^2)-v*atan((y-mu_pearson)/sigma_pearson));
        end
        %%%% Taper Prior
        if ismember(param.Method,[10,11,12,13])
            Taper_return = Taper(b);
            C = C + Taper_return{1};
            if (param.Use_Gradient == 1)
                grad_arr = grad_arr + Taper_return{2};
            end
        end
        %%%% Gaussian/Truncated Gaussian Prior
        if (~isempty(param.PriorCell) || ~isempty(Prior_Current) || (param.Method == 13))
            if ~isempty(Prior_Current) 
                Prior_LT = Prior_Current{1};
            else
                Prior_LT = param.PriorCell{1};
            end
            
            for idx_tau2 = 1:param.order
                Prior_LT_current = Prior_LT{idx_tau2};
                if (ismember(idx_tau2,param.Prior_order))
                    if (param.randPrior)
                        Prior_LT_current = Prior_LT_current + rand_arr(idx_tau2)*param.sigma_prior; 
                    end

                    if (param.Prior_truncated == 0)
                        Reg = param.RegFactor1*((Get_param(b,param,"tau",idx_tau2) - Prior_LT_current).^2) / (param.sigma_prior^2);
                        C = C + Reg;
                        if (param.Use_Gradient == 1)
                            Reg_grad = 2*param.RegFactor1*((Get_param(b,param,"tau",idx_tau2) - Prior_LT_current)) / (2 * param.sigma_prior^2);
                            grad_arr(Get_param(b,param,"tau",-idx_tau2)) = grad_arr(Get_param(b,param,"tau",-idx_tau2)) + Reg_grad;
                        end
                    else
                        t_mean = Prior_LT_current;
                        t_std = param.sigma_prior;
                        tau_l = t_mean-t_std;
                        tau_u = t_mean+t_std;
                        tau_vec = param.tau_vec;
                        rightIdx = tau_vec > tau_u;
                        middleIdx = (tau_vec >= tau_l) & (tau_vec <= tau_u);
                        leftIdx = tau_vec < tau_l;
                        Reg = zeros(size(param.tau_vec)); % initialize
                        Reg(leftIdx) = -((tau_vec(leftIdx) - tau_l).^2) / (2 * t_std^2);
                        Reg(middleIdx) = 0;                
                        Reg(rightIdx) = -((tau_vec(rightIdx) - tau_u).^2) / (2 * t_std^2);
                        Reg = 5000*Reg;
                        [~,tau_idx] = min(abs(param.tau_vec - Get_param(b,param,"tau",idx_tau2)));
                        if (param.Use_Gradient == 1)
                            Reg_grad = zeros(length(Reg),1);
                            Reg_grad(leftIdx) = (tau_vec(leftIdx) - tau_l)/(2*t_std^2);
                            Reg_grad(middleIdx) = 0;
                            Reg_grad(rightIdx) = (tau_vec(rightIdx) - tau_u)/(2*t_std^2);
                            grad_arr(Get_param(b,param,"tau",-idx_tau2)) = grad_arr(:,Get_param(b,param,"tau",-idx_tau2)) + Reg_grad(tau_idx);
                        end
                        C = C + Reg(tau_idx);
                    end
                end
            end
        end

        if (param.fixed_amp_prior == 1)
            for idx_tau3 = 1:param.order
                if (param.AmpEst_flag == 1)
                    Prior_temp  = Prior_Current{2};
                    Prior_Amp = Prior_temp{idx_tau3};
                else
                    Prior_Amp = param.PriorCell{2}{idx_tau3};
                end
                Prior_Amp_current = Prior_Amp;
                if (param.randPrior)
                    Prior_Amp_current = Prior_Amp + rand_arr(param.order+idx_tau3)*param.sigma_prior_amp; 
                end
                Reg = param.RegFactor_amp*((Get_param(b,param,"amp",idx_tau3) - Prior_Amp_current).^2) / (param.sigma_prior_amp^2);
                C = C + Reg;
                if (param.Use_Gradient == 1)
                    Reg_grad = 2*param.RegFactor_amp*((Get_param(b,param,"amp",idx_tau3) - Prior_Amp_current)) / (2 * param.sigma_prior_amp^2);
                    grad_arr(Get_param(b,param,"tau",-idx_tau3)) = grad_arr(Get_param(b,param,"tau",-idx_tau3)) + Reg_grad;
                end
            end
        end


        function Taper_return = Taper(B)
            Cost_taper = 0;
            %%% Lifetime taper
            if (param.Use_Gradient == 1)
                taper_grad = zeros(1,length(B));
            end
            tau_axis = param.tau_vec;
            taper_axis = zeros(length(tau_axis),1);
            tau_lower_subaxis_idx = tau_axis < param.taper_tau_lower;
            tau_upper_subaxis_idx = tau_axis > param.taper_tau_upper;
            taper_axis(tau_lower_subaxis_idx) = param.taper_regl*(tau_axis(tau_lower_subaxis_idx) - param.taper_tl).^2;
            taper_axis(tau_upper_subaxis_idx) = param.taper_regu*(tau_axis(tau_upper_subaxis_idx) - param.taper_tu).^2;

            if (param.Use_Gradient == 1)
                taper_grad_axis = zeros(length(tau_axis),1);
                taper_grad_axis(tau_lower_subaxis_idx) = 2*param.taper_regl*(tau_axis(tau_lower_subaxis_idx) - param.taper_tl);
                taper_grad_axis(tau_upper_subaxis_idx) = 2*param.taper_regu*(tau_axis(tau_upper_subaxis_idx) - param.taper_tu);
            end

            for idx_tau4 = 1:param.order
                [~,tau_axis_idx] = min(abs(tau_axis - Get_param(B,param,"tau",idx_tau4)));
                if (param.Use_Gradient == 1)
                    taper_grad(Get_param(B,param,"tau",-idx_tau4)) = taper_grad_axis(tau_axis_idx);
                end
                Cost_taper = Cost_taper + taper_axis(tau_axis_idx);
            end

            %%% Amplitude taper
            if (param.order > 1 && param.fixed_amp == 0)
                amp_axis = 0:0.01:1;
                amp_subaxis_idx_lower = amp_axis < param.taper_amp_limit;
                amp_subaxis_idx_upper = amp_axis > (1-param.taper_amp_limit);

                taper_axis = zeros(length(amp_axis),1);               
                taper_axis(amp_subaxis_idx_lower) = param.Ataper_reg*(amp_axis(amp_subaxis_idx_lower) - param.taper_amp_limit).^2;
                taper_axis(amp_subaxis_idx_upper) = param.Ataper_reg*(amp_axis(amp_subaxis_idx_upper) - (1-param.taper_amp_limit)).^2;

                if (param.Use_Gradient == 1)
                    taper_grad_axis(amp_subaxis_idx_lower) = 2*param.Ataper_reg*(amp_axis(amp_subaxis_idx_lower) - param.taper_amp_limit);
                    taper_grad_axis(amp_subaxis_idx_upper) = 2*param.Ataper_reg*(amp_axis(amp_subaxis_idx_upper) - (1-param.taper_amp_limit));                
                end
                
                for idx_tau5 = 1:param.order
                    [~,amp_axis_idx] = min(abs(amp_axis - Get_param(B,param,"amp",1)));
                    if (param.Use_Gradient == 1)
                        taper_grad(Get_param(B,param,"amp",-idx_tau5)) = taper_grad_axis(amp_axis_idx);
                    end
                    Cost_taper = Cost_taper + taper_axis(amp_axis_idx);
                end
            end
            Taper_return{1} = Cost_taper;
            if param.Use_Gradient == 1
                Taper_return{2} = taper_grad;
            end
        end
        varargout{1} = C;
        if (param.Use_Gradient == 1)
            % grad_arr = grad_arr / norm(grad_arr);
            % disp(grad_arr);
            varargout{2} = grad_arr;
        end
    end

    function plot_I (~,~) 
        D = (y_fit - y_OG)./sqrt(y_fit);
        exp_dist = param.tau_dist;
        
        gcf = figure(7);sgtitle(strcat('Lifetime Fitting, \chi^{2}_{r} = ', num2str(chi_r),' ->', num2str(chi_r_clipped)));
        set(gcf,'position',[100,-10,700,700])
        ax1 = subplot(5,1,[1 2]);hold on;
        yyaxis left
            plot(t,y_OG,'b','linewidth',2);
            plot(t,y_fit,'k-.','linewidth',1.5);
            set(gca, 'YScale', 'log')
        grid on;
        ylabel('Intensity','fontweight', 'bold','fontsize', fontsz);            
        
        yyaxis right
            plot(t_IRF,H,'r','linewidth',1.5);

        legend('Measured','Fitted');
        ax2 = subplot(5,1,3);hold on;
        plot(t,D,'k','linewidth',2)
        ylabel('Deviation','fontweight', 'bold','fontsize', fontsz);            
        grid on;
        % ylim([-10,10]);

        ax4 = subplot(5,1,4);hold on; 
        plot(t,chi_vec,'k','linewidth',2,'Linestyle','none','Marker','.')
        ylabel({"\chi^{2}";"Contribution"},'fontweight', 'bold','fontsize', fontsz);            
        grid on;
        xlabel('Time (ns)','fontweight', 'bold','fontsize', fontsz);
        linkaxes([ax1,ax2,ax4],'x')
        % xlim([t(idx_peak+3),max(t)])
        xlim([0,max(t)])
        subplot(5,1,5);hold on; 
        if (param.cont == 0)
            exp_dist (exp_dist == 0 ) = NaN;
            stem(param.tau_vec,exp_dist,'k','linewidth',2)
        else
            plot(param.tau_vec,param.tau_dist,'k','linewidth',2);
        end
        grid on;
        xlim([0,param.tau_vec(end)]);
        if (param.cont == 0)
            Str_type = 'Discrete';
        elseif (param.cont_type == 1)
            Str_type = 'Gaussian';
        elseif (param.cont_type == 2)
            Str_type = 'Gamma';
        elseif (param.cont_type == 3)
            Str_type = 'Inverse Gamma';
        elseif (param.cont_type == 4)
            Str_type = 'Log-normal';
        end
        ylabel({"Amplitude";Str_type},'fontweight', 'bold','fontsize', fontsz);            
        xlabel('Lifetime (ns)','fontweight', 'bold','fontsize', fontsz);
        

        if(param.PLOT_Tau)   %%% Plot frequency info
            figure(72);hold on
            if (param.cont == 0)
                exp_dist (exp_dist == 0 ) = NaN;
                stem(param.tau_vec,exp_dist,'linewidth',2)
            else
                plot(param.tau_vec,param.tau_dist,'linewidth',2);
            end
            grid on;
            xlim([0,param.tau_vec(end)]);
            ylabel('Amplitude','fontweight', 'bold','fontsize', fontsz);            
            xlabel('Lifetime (ns)','fontweight', 'bold','fontsize', fontsz);
        end
        if(param.PLOT_Freq)   %%% Plot frequency info
            y_W = FT(y,dt_IRF);
            H_W = FT(H,dt_IRF);
            output_w = @(b) FT(output_t(b),dt_IRF);
            figure(8);hold on;
            plot(f,abs(y_W)/max(abs(y_W)),'b','linewidth',2);
            plot(f,abs(H_W)/max(abs(H_W)),'r-.','linewidth',2);
            plot(f,abs(output_w(B))/max(abs(output_w(B))),'k-.','linewidth',2);
            grid on;
            xlabel('Freq (GHz)');
            ylabel('Intensity');
            legend('Measured','IRF','Fitted');
        end
        % fig = figure(17);clf;hold on;grid on
        % if (param.NeighborsNoise_flag == 1 || param.Method == 7)
        %     plot(t,NoiseVector_clipped,'b','LineWidth',2,"Displayname","Neighbors Noise")
        % end
        % plot(t,NoiseVector_empirical_clipped,'Color',[0.47,0.67,0.19],'LineWidth',2,"Displayname","S-Empirical Model")
        % plot(t,NoiseVector_empirical_neighbors_clipped,'Color',[1.00,0.41,0.16],'LineWidth',2,"Displayname","M-Empirical Model")
        % % plot(t,NoiseActual_clipped,'k','LineWidth',2,"Displayname","Simulated Noise")
        % plot(t,y_OG,'Color',[0.00,0.45,0.74],'LineWidth',2,"Displayname","Poisson Model")
        % % title("Signal and Noise")
        % xlim([t(1),t(end)])
        % ylim([10,2e5])
        % xlabel("Time (ns)")
        % ylabel("Noise Variance")
        % set(gca, 'YScale', 'log')
        % axx = gca;
        % axx.FontSize = 14;
        % pos = fig.Position;
        % fig.Position = [pos(1) pos(2) 400 350];
        % legend

    end

    function exp_dist = get_dist(B)
        tau_vector = param.tau_vec;
        exp_dist = zeros(size(tau_vector));
        if (param.SolveMatrix == 1)
            [~,A] = output_t(B);
        end
        if (param.cont == 0)
            for j = 1:param.order
                idx = find(abs(tau_vector - Get_param(B,param,"tau",j)) == min(abs(tau_vector - Get_param(B,param,"tau",j))));
                exp_dist(idx) = exp_dist(idx) + Get_param(B,param,"amp",j);
            end
        else
            for j = 1:param.order
                mu = Get_param(B,param,"tau",j);
                sigma = Get_param(B,param,"sigma",j);
                if(param.cont_type == 1)    %%%% Gaussian distribution
                    pdf_temp = normpdf(tau_vector,mu,sigma)*(tau_vector(end)-tau_vector(1))/length(tau_vector);
                elseif(param.cont_type == 2)    %%%% Gama distribution
                    pdf_temp = gampdf(tau_vector,mu^2/sigma,sigma/mu)*(tau_vector(end)-tau_vector(1))/length(tau_vector);
                elseif(param.cont_type == 3)    %%%% Inverse Gama distribution
                    pdf_temp = gampdf(1./tau_vector, mu^2 / sigma, sigma/(mu^3)) * (tau_vector(end) - tau_vector(1)) / length(tau_vector);    
                elseif(param.cont_type == 4)    %%%% Log-normal distribution
                    fn_gamma = asinh(sigma/(2*mu));
                    pdf_temp = exp(-((log(1./tau_vector)-log(1/mu))/fn_gamma).^2);    
                end
                exp_dist = exp_dist + Get_param(B,param,"amp",j)*pdf_temp;
            end
        end
    end
    
    function [problem,IC] = get_problem(~,~)
        if (param.SolveMatrix == 0)
            lb = [];
            ub = [];
            IC = [];
            Aeq = [];
            Beq = 1;
            Aineq = []; 
            if (param.order == 1)
                Bineq = [];
            else
                Bineq = 0.01*ones(param.order-1,1);
            end
            if (param.tshift_flag == 1)
                lb = [lb, -param.t_shift_max];
                ub = [ub, param.t_shift_max];
                IC = [IC, 0];
                Aeq = [Aeq, 0];
                if (param.order > 1)
                    Aineq = [Aineq,zeros(param.order-1,1)];
                end
            end
            if (param.Scaleshift_flag == 1)
                lb = [lb, 1];
                ub = [ub, inf];
                IC = [IC, Ar_IC];
                Aeq = [Aeq, 0];
                if (param.order > 1)
                    Aineq = [Aineq,zeros(param.order-1,1)];
                end
            end

            if (param.cont == 0)
            %%% Multi-exponential
            %%% Parameters are: [Shift  Amp   tau-1  Amp-1 ... ]
            %%% Parameters are: [b(1)   b(2)  b(3)   b(4)  ... ]
                for k = 1:param.order
                    if (param.fixed_amp == 1 || param.order == 1)
                        lb = [lb,param.lb];
                        ub = [ub,param.ub];
                        IC = [IC,k*param.IC];
                        if (param.order > 1)
                            Aineq = [Aineq,zeros(param.order-1,1)];
                        end
                    else
                        lb = [lb,param.lb,param.amplb];
                        ub = [ub,param.ub,param.ampub];
                        IC = [IC,k*param.IC,1/param.order];
                        Aeq =[Aeq,0,1];
                        if (param.order > 1)
                            Aineq = [Aineq,zeros(param.order-1,2)];
                        end
                    end
                    if (param.order > 1)
                        if(k ~= 1)
                            Aineq(k-1,Get_param(IC, param, "tau", -k)) = -1;
                            Aineq(k-1,Get_param(IC, param, "tau", -(k-1))) = 1;
                        end
                    end
                end
            elseif (param.cont == 1)
                for k = 1:param.order
                    if (param.fixed_amp == 1 || param.order == 1)
                        lb = [lb,param.lb,param.sigma_max];
                        ub = [ub,param.ub,param.sigma_max];
                        IC = [IC,k*param.IC,param.sigma_max];
                    else
                        lb = [lb,param.lb,param.sigma_max,param.amplb];
                        ub = [ub,param.ub,param.sigma_max,param.ampub];
                        IC = [IC,k*param.IC,param.sigma_max,1/param.order];
                        Aeq =[Aeq,0,0,1];
                    end
                end
            end
            % %%% Range of the added DC is set here
            if (param.DCshift_flag == 1)
                lb = [lb,-5];
                ub = [ub,5];
                IC = [IC,0];
                Aeq =[Aeq,0];
                if (param.order > 1)
                    Aineq = [Aineq,zeros(param.order-1,1)];
                end
            end

        else
            lb = [];
            ub = [];
            IC = [];
            Aeq =[];

            if (param.cont == 0)
                for k = 1:param.order
                    lb = [lb,param.lb];
                    ub = [ub,param.ub];
                    IC = [IC,k*param.IC];
                end
                if (param.tshift_flag == 1)
                    lb = [lb,-param.t_shift_max];
                    ub = [ub,param.t_shift_max];
                    IC = [IC,0];
                end
                %%%
                if (param.Method == 5)
                    lb = [lb,0];
                    ub = [ub,10];
                    IC = [IC,1];
                end
                if (param.DCshift_flag == 1)
                    lb = [lb,-DarkCurrent];
                    ub = [ub,DarkCurrent];
                    IC = [IC,0];
                end
            end
        end
        opts = optimoptions('fmincon', 'Algorithm','interior-point','MaxFunEvals', 50000, 'MaxIter', 5000, 'StepTolerance', 1e-10,'Display', 'off');

        if (param.Use_Gradient == 1)
            opts.SpecifyObjectiveGradient = true;
        end
        if (param.fixed_amp == 0 && param.order > 1)
            problem = createOptimProblem('fmincon','x0',IC,'objective',@Costfn,'lb',lb,'ub',ub,'Aineq',Aineq,'bineq',Bineq,'Aeq', Aeq, 'beq', Beq,'options', opts);
        elseif (param.fixed_amp == 1 || param.order == 1)
            problem = createOptimProblem('fmincon','x0',IC,'objective',@Costfn,'lb',lb,'ub',ub,'Aineq',Aineq,'bineq',Bineq,'Aeq', [], 'beq', [],'options', opts);
        end
    end
     
    function varargout = output_t(B)
        %%% Parameters are: [shift Amp  mu   sigma amp  ... ]
        %%% Parameters are: [b(1)  b(2) b(3) b(4)  b(5) ...]
        if (param.SolveMatrix == 0)
            exp_vector = zeros(size(t_IRF,1),1);
            tau_vector  = [cast(1000*Get_param(B,param,"tau",1),"int32")];
            if (param.cont == 1)
                for k = 2:param.order                    
                    if (ismember(cast(1000*Get_param(B,param,"tau",k),"int32"),tau_vector))
                        B(Get_param(B,param,"tau",-k)) = Get_param(B,param,"tau",k)+0.001;
                        tau_vector  = [tau_vector cast(1000*Get_param(B,param,"tau",k),"int32")];
                    else
                        tau_vector  = [tau_vector cast(1000*Get_param(B,param,"tau",k),"int32")];
                    end
                end
                exp_dist = get_dist(B);
                exp_vector = sum(exp_dist .* exp(-t_IRF./ param.tau_vec),2);
            end
            if (param.cont == 0)
                for k = 2:param.order
                    if (ismember(cast(1000*Get_param(B,param,"tau",k),"int32"),tau_vector))
                        B(Get_param(B,param,"tau",-k)) = Get_param(B,param,"tau",k)+0.001;
                    end
                    tau_vector  = [tau_vector cast(1000*Get_param(B,param,"tau",k),"int32")];
                end
                for k = 1:param.order
                    tau_temp = Get_param(B,param,"tau",k);
                    if (param.NormDivTau == 1)
                        exp_vector = exp_vector + Get_param(B,param,"amp",k)*exp(-t_IRF/tau_temp)/tau_temp;
                    else
                        exp_vector = exp_vector + Get_param(B,param,"amp",k)*exp(-t_IRF/tau_temp);
                    end
                end
            end
            exp_vector(isnan(exp_vector))=0;    %%% To avoid any errors when first point is NaN at 0 lifetime            
            if (param.tshift_flag == 1)
                output_temp  = myconv(H,exp_vector,Get_param(B,param,"tshift"));
            else
                output_temp = myconv(H,exp_vector,0);
            end
            output_temp2 = Get_param(B,param,"Ar")*output_temp;
            output  = abs(output_temp2) + DarkCurrent;
            if (param.DCshift_flag == 1)
                output  = output + Get_param(B,param,"DC");
            end
            % output(output>4095) = 4095;
            output = output(common_time);
            varargout{1} = output(clip_vector);
            %%% Gradients
            if (param.Use_Gradient == 1)
                grad_arr = zeros(length(output(common_time)),length(B));
                if (param.tshift_flag == 1) 
                    output_W = FT(cat(1,output_temp,zeros(N_ZP-N,1)), dt_IRF);
                    dW_dtshift = -1i * 2 * pi * f_ZP .* output_W;
                    d_output_dtshift = real(ifftshift(IFT(dW_dtshift, 1/(N_ZP*dt_IRF))));
                    d_output_dtshift = Ar_const*Get_param(B,param,"Ar")*d_output_dtshift(1:length(output_temp));
                    grad_arr(:,Get_param(B,param,"tshift",-1)) = d_output_dtshift(common_time);
                end
                grad_arr(:,Get_param(B,param,"Ar",-1)) = output_temp(common_time);
                if (param.DCshift_flag == 1) 
                    grad_arr(:,Get_param(B,param,"DC",-1)) = ones(1,length(output_temp(common_time)));
                end
                for idx_order_2 = 1:param.order
                    tau_temp = Get_param(B,param,"tau",idx_order_2);
                    if (param.NormDivTau == 0)
                        output_temp_tau = myconv(H,((t_IRF)/tau_temp^2).*exp(-t_IRF/tau_temp),Get_param(B,param,"tshift"));
                        grad_arr(:,Get_param(B,param,"tau",-idx_order_2)) = Get_param(B,param,"amp",idx_order_2)*Ar_const*Get_param(B,param,"Ar")*output_temp_tau(common_time);
                    else %%% This needs to be fixed
                        output_temp_tau = myconv(H,((t_IRF-tau_temp)/tau_temp^3).*exp(-t_IRF/tau_temp),Get_param(B,param,"tshift"));
                        grad_arr(:,Get_param(B,param,"tau",-idx_order_2)) = Get_param(B,param,"amp",idx_order_2)*Ar_const*Get_param(B,param,"Ar")*output_temp_tau(common_time);

                    end
                    
                    if (param.order > 1 && param.fixed_amp == 0)
                        if (param.NormDivTau == 0)
                            output_temp_amp = myconv(H,exp(-t_IRF/tau_temp),Get_param(B,param,"tshift"));
                        else 
                            output_temp_amp = myconv(H,exp(-t_IRF/tau_temp)/tau_temp,Get_param(B,param,"tshift"));
                        end
                        grad_arr(:,Get_param(B,param,"amp",-idx_order_2)) = Ar_const*Get_param(B,param,"Ar")*output_temp_amp(common_time);
                    end
                end
                varargout{2} = grad_arr;
            end
        end
        if (param.SolveMatrix == 1)
            tau_vector  = [cast(1000*Get_param(B,param,"tau",1),"int32")];
            for k = 2:param.order
                if (ismember(cast(1000*Get_param(B,param,"tau",k),"int32"),tau_vector))
                    B(Get_param(B,param,"tau",-k)) = Get_param(B,param,"tau",k)+0.001;
                    tau_vector  = [tau_vector cast(1000*Get_param(B,param,"tau",k),"int32")];
                else
                    tau_vector  = [tau_vector cast(1000*Get_param(B,param,"tau",k),"int32")];
                end
            end
            for k = 1:param.order
                if (param.tshift_flag == 1)
                    exp_temp = myconv(H,exp(-t_IRF/Get_param(B,param,"tau",k)),Get_param(B,param,"tshift"))+DarkCurrent+Get_param(B,param,"DC");
                else
                    exp_temp = myconv(H,exp(-t_IRF/Get_param(B,param,"tau",k)),0)+DarkCurrent+Get_param(B,param,"DC");
                end
                exp_temp = exp_temp(common_time);
                exp_trim{k} = exp_temp(clip_vector);
                exp_full{k} = exp_temp;
            end            

            if (param.Method == 6)
                 variance_NoTrim = NoiseVector_empirical;
                 variance_Trim   = NoiseVector_empirical_clipped;
            elseif (param.Method == 7)
                 variance_NoTrim = NoiseVector;
                 variance_Trim   = NoiseVector_clipped;
            elseif (param.Method == 8)
                 variance_NoTrim = NoiseActual;
                 variance_Trim   = NoiseActual_clipped;
            elseif(param.Method == 9)
                 variance_NoTrim = NoiseVector_empirical_neighbors;
                 variance_Trim   = NoiseVector_empirical_neighbors_clipped;
            elseif(param.Method == 2)
                 variance_NoTrim = param.MCPNoise(y_OG);
                 variance_Trim   = variance_NoTrim(clip_vector);
            else
                 variance_NoTrim = y_OG;
                 variance_Trim = y_OG(clip_vector);
            end
            covm = diag(variance_Trim);
            fitmat = [ones(size(variance_Trim))];
            outmat = [ones(size(variance_NoTrim))];
            for k2 = 1:param.order
                fitmat = [fitmat ,reshape(exp_trim{k2},size(variance_Trim))];    
                outmat = [outmat ,reshape(exp_full{k2},size(variance_NoTrim))];    
            end
            A = lscov(fitmat,reshape(y(clip_vector),size(variance_Trim)),covm);
            varargout{1} = outmat*A;
            varargout{2} = A;
        end
    end

    function C = myconvshift(A,t_shift,dt)   % perform 1D Linear discrete convolution followed by clipping
        N = length(A);
        f_max =1/dt;
        df=1/(N*dt);
        f = (-f_max/2:df:f_max/2-df)';
        A_W = FT(A,dt);
        C_W = A_W.*exp(-1j*2*pi*f*t_shift);
        C = abs(IFT(C_W,1/(N*dt)));
    end

    function C = mytimeshift(A, t_shift, dt)
        shift_amount = round(t_shift / dt);
        C = circshift(A, shift_amount);
    end

    function C = myconv(A, B,t_shift)   % perform 1D Linear discrete convolution followed by clipping
        N = length(A);
        %%% Zero Padding
        A_ZP = cat(1,A,zeros(N_ZP-N,1));
        B_ZP = cat(1,B,zeros(N_ZP-N,1));
        %%% FFT --> shift and Conv --> IFFT
        A_W = FT(A_ZP,dt_IRF);
        B_W = FT(B_ZP,dt_IRF);
        try
            C_W = A_W.*B_W.*exp(-1j*2*pi*f_ZP*t_shift);
            % C_try = conv(A_ZP, B_ZP);
        catch ME
            err_str = ME;
            if (param.DEBUGGING == 1)
                keyboard;
            else
                err_status = 1;
                if (length(A_W)==length(B_W))
                    err_status = 2;
                elseif (length(f_ZP)==length(A_W))
                    err_status = 3;
                elseif (length(f_ZP)==length(B_W))
                    err_status = 4;
                end
                return;
           end
        end
        %%% Salem, this might need to be N_zp. And the scaling might be
        %%% removed all together
        C = real(IFT(C_W,1/(N_ZP*dt_IRF)));
        %%% Clipping the zero padded signal
        C = C(1:N);
    end
    
    function G = FT(g, dt)  % 1D DFT where delta is sample width
        G = fftshift(fft(fftshift(g))) * dt;
    end
    
    function g = IFT(G, df) % 1D DIFT where df is sample with in frequency domain (1/N/dt)
        g = (ifft(ifftshift(G)) * length(G) * df);
    
    end

    function Debug_plots(~,~) 
        figure(1);hold on;
        plot(t_IRF,myconv(H,exp(-t_IRF/2),-1  ),'b');
        plot(t_IRF,myconv(H,exp(-t_IRF/2),-0.9),'b');
        plot(t_IRF,myconv(H,exp(-t_IRF/2),-0.8),'b');
        plot(t_IRF,myconv(H,exp(-t_IRF/2),-0.7),'b');
        plot(t_IRF,myconv(H,exp(-t_IRF/2),-0.6),'b');
        plot(t_IRF,myconv(H,exp(-t_IRF/2),-0.5),'b');
        plot(t_IRF,myconv(H,exp(-t_IRF/2),-0.4),'b');
        plot(t_IRF,myconv(H,exp(-t_IRF/2),-0.3),'b');
        plot(t_IRF,myconv(H,exp(-t_IRF/2),-0.2),'b');
        plot(t_IRF,myconv(H,exp(-t_IRF/2),-0.1),'b');
        plot(t_IRF,myconv(H,exp(-t_IRF/2), 0  ),'k','linewidth',2);
        plot(t_IRF,myconv(H,exp(-t_IRF/2), 0.1),'r');
        plot(t_IRF,myconv(H,exp(-t_IRF/2), 0.2),'r');
        plot(t_IRF,myconv(H,exp(-t_IRF/2), 0.3),'r');
        plot(t_IRF,myconv(H,exp(-t_IRF/2), 0.4),'r');
        plot(t_IRF,myconv(H,exp(-t_IRF/2), 0.5),'r');
        plot(t_IRF,myconv(H,exp(-t_IRF/2), 0.6),'r'); 
        plot(t_IRF,myconv(H,exp(-t_IRF/2), 1  ),'r');
        keyboard 
    end

      function [Ar_IC,tshift] = align_rising_edge(meas, temp, DarkCurrent)
        [~, peak_idx_meas] = max(meas);
        [~, peak_idx_temp] = max(temp);
        tshift = dt_sig*peak_idx_meas - dt_IRF*peak_idx_temp;
        Ar_IC = (max(meas)-DarkCurrent)/max(temp);
    end
    
    function B_selected = Get_param(B,param,param_name,idx)
        
        %%% Discrete Parameters are: [tshift A_r tau1 Amp1  tau2 amp2 ... DC]
        %%% Continious Parameters are: [tshift A_r mu1 sigma1 amp1 ... DC]
    
        if (strcmp(param_name,"DC"))
            if (param.DCshift_flag == 1)
                B_selected = B(end);
                if (nargin > 3 && idx == -1)
                    B_selected = length(B);
                end
            else 
                B_selected = 0;
            end
            return;
        elseif (strcmp(param_name,"tshift"))
            if (param.tshift_flag == 1)
                B_selected = B(1) + tshift;
                if (nargin > 3 && idx == -1)
                    B_selected = 1;
                end
            else 
                B_selected = 0;
            end
            return;
        elseif (strcmp(param_name,"Ar"))
            if (param.Scaleshift_flag == 1)
                if (param.tshift_flag == 1)
                    B_selected = B(2);
                    if (nargin > 3 && idx == -1)
                        B_selected = 2;
                    end
                else 
                    B_selected = B(1);
                    if (nargin > 3 && idx == -1)
                        B_selected = 1;
                    end
                end
            else
                B_selected = 1;
            end
            return;
        elseif (strcmp(param_name,"amp") && param.order == 1)
            B_selected = 1;
            return;
        end
            
        const_param_count = 0; 
        if (param.tshift_flag)
            const_param_count = const_param_count + 1; 
        end
        if (param.Scaleshift_flag)
            const_param_count = const_param_count + 1; 
        end
        if (param.fixed_amp == 0)
            tau_idx = zeros(param.order,1);
            for idx_order = 1:param.order
                if (param.cont == 0)
                    tau_idx(idx_order) = const_param_count+2*(idx_order-1)+1;
                else
                    tau_idx(idx_order) = const_param_count+3*(idx_order-1)+1;
                end
            end
            if (param.cont == 1)
                sigma_values = B(tau_idx+1);
                amp_idx = tau_idx+2;
            else
                amp_idx = tau_idx+1;
            end
            if (param.order > 1)
                amp_value = B(amp_idx(abs(idx)));
            end
            tau_value = B(tau_idx(abs(idx)));
        else
            tau_idx = zeros(param.order,1);
            for idx_order = 1:param.order
                if (param.cont == 0)
                    tau_idx(idx_order) = const_param_count+1+(idx_order-1);
                else
                    tau_idx(idx_order) = const_param_count+1+2*(idx_order-1);
                end
            end
            tau_value = B(tau_idx(abs(idx)));
            amp_value = param.fixed_amp_arr;
            if (param.cont == 1)
                sigma_values = B(tau_idx+1);
            end
        end
    
        switch param_name 
        case "tau"
            if (idx > 0)
                B_selected = tau_value;
            elseif (idx < 0)
                B_selected = tau_idx(abs(idx));
            end
            return;
        case "amp"
            if (param.fixed_amp == 0)
                if (idx > 0)
                    B_selected = amp_value;
                elseif (idx < 0)
                    B_selected = tau_idx(abs(idx))+1;
                end
            else
                B_selected = amp_value(abs(idx));
            end
            return;
    
        case "sigma"
            B_selected = sigma_values(idx);
            if (B_selected == 0)
                B_selected = 0.001;
            end
            return;
        end
    end

function [tau] = phasor_analysis(y_raw, t,f,tau_ref)
    tau_ref = tau_ref*1e-9;
    y_temp = y_raw-DarkCurrent;
    omega = 2 * pi * f;
    % y(y < 0) = 0;  % Clip negative values
    
    [M,idx_IRFmax] = max(H(common_time));
    y_temp = y_temp(idx_IRFmax:end);
    y_temp = y_temp / sum(y_temp);  % area normalization to avoid scaling effects
    t_use = t(idx_IRFmax:end)*1e-9;
    t_use = t_use-t_use(1);

    y_ref = y_Actual-DarkCurrent;
    % y_ref(y_ref < 0) = 0;  % Clip negative values
    y_ref = y_ref(idx_IRFmax:end);
    y_ref = y_ref / sum(y_ref);  % area normalization to avoid scaling effects
    g_ref = sum(y_ref.*cos(omega*t_use));
    s_ref = sum(y_ref.*sin(omega*t_use));
    g_theo = 1 / (1 + (omega * tau_ref)^2);
    s_theo = (omega * tau_ref) / (1 + (omega * tau_ref)^2);
    phi_meas = atan2(s_ref, g_ref);         % measured phase
    phi_theo = atan2(s_theo, g_theo);       % theoretical phase
    delta_phi = phi_theo - phi_meas;        % additive phase correction (in radians)

    M_meas = sqrt(g_ref^2 + s_ref^2);       % measured modulation
    M_theo = sqrt(g_theo^2 + s_theo^2);     % theoretical modulation
    Mscale = M_theo / M_meas;               % multiplicative modulation correction

    
    g = sum(y_temp.*cos(omega*t_use));
    s = sum(y_temp.*sin(omega*t_use));
    g_corr = Mscale * (g * cos(delta_phi) - s * sin(delta_phi));
    s_corr = Mscale * (g * sin(delta_phi) + s * cos(delta_phi));

    tau = 1e9*s_corr/(g_corr*omega);
end



end 



function y_avg = FindDiff(y, N)
    L = length(y);
    y_avg = zeros(size(y)); % preallocate
    y_avg(1) = y(1);
    y_avg(L) = y(L);
    for i = 2:L-1
        idx_start = max(1, i - N);
        idx_end   = min(L, i + N);
        neighbor_indices = idx_start:idx_end;
        neighbor_indices(neighbor_indices == i) = []; % Exclude current point
        y_avg(i) = mean(y(neighbor_indices));
    end
end

function FixFilter(param, path_files_temp, filename_exist_arr, Mask_OG, Skip_OG)
    N = numel(path_files_temp);

    % Helper to read peak frame
    function Ipeak = read_peak_frame(pth, param)
        pth = char(pth);
        files = dir([pth, '/*.tif']);
        if isempty(files)
            Ipeak = [];
            return;
        end
        list_raw = split(strip(ls([pth, '/*.tif'])));
        list_raw = list_raw(~cellfun('isempty',list_raw));
        if ~isempty(param) && isfield(param,'idx_peak_thresh') && ~isempty(param.idx_peak_thresh)
            idxp = param.idx_peak_thresh;
            fn = erase(string(list_raw{idxp}),"'");
            th = Tiff(fn,'r');
            Ipeak = double(read(th));
        else
            % auto-peak by sum over pixels
            paramsFile = fullfile(pth,'RecSettings.txt');
            [~,delta_t,tmin,tmax,~,BINNING,~,~] = GetExperimentalParameters(paramsFile); %#ok<ASGLU>
            t_sig = (tmin:delta_t:tmax)'; %#ok<NASGU>
            n = length(t_sig);
            I_temp = zeros(1216/BINNING, 1936/BINNING, n);
            for ii = 1:n
                fn = erase(string(list_raw{ii}),"'");
                th = Tiff(fn,'r');
                I_temp(:,:,ii) = double(read(th));
            end
            Isum = squeeze(sum(I_temp,[1 2]));
            [~,idxp] = max(Isum);
            Ipeak = I_temp(:,:,idxp);
        end
    end
    satMask = cell(N,1);
    Mask_out  = cell(N,1);
    Mask_all = [];
    for k = 1:N
        if ~filename_exist_arr(k), continue; end
        if (Skip_OG{k}), continue; end
        pth = path_files_temp{k};
        Ipeak = read_peak_frame(pth, param);    % Load peak frame
        satMask{k} = Ipeak >= 4095; 
        Mask_out{k} =  Mask_OG{k}.*~satMask{k};
        if isempty(Mask_all)
            Mask_all = Mask_out{k};            
        else
            Mask_all = Mask_all.*Mask_out{k};
        end 
    end
    for k = 1:N
        if Skip_OG{k}, continue; end
        if ~filename_exist_arr(k), continue; end
        if (param.Mask_together)
            Mask = Mask_all;
        else
            Mask = Mask_out{k};
        end
        Mask_filepath = Get_savefile_name(param,path_files_temp{k},"Mask");
        save(Mask_filepath,'Mask','-append');
    end
end


