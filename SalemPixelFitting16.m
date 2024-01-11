function status = SalemPixelFitting16(arg1, param,param_plot)

    path = char(arg1);
    status = 1;
    %%% Setting file name
    if (path(end) == '/')
        path = path(1:end-1);
    end
    filename_saving_temp = split(strip((path)),'/');
    filename_saving = strcat(string(filename_saving_temp(end-1)),string(filename_saving_temp(end)));
    if (param.cont == 1)
        str_cont = 'C';
    else
        str_cont = 'D';
    end
    filename_saving_name = strcat('Results/RESULTS_EXP',num2str(param.order),str_cont,'_',filename_saving,'_M',num2str(param.Method));
    filename_saving_results = strcat(filename_saving_name,'.mat');
    %%% Extracting file parameters
    if (param.SHOW_ONLY == 0)
        parametersFile = strcat(path,'/RecSettings.txt'); % contains experimental parameters
        [~,delta_t,tmin,tmax,Gate,BINNING,MCP] = GetExperimentalParameters(parametersFile);    % get experimental parameters
        param.mcp = MCP;
        HorizontalPixels = 1:1936/BINNING;     % based on binning factor [pixels]
        VerticalPixels = 1:1216/BINNING;         % based on binning factor [pixels]
        t_sig = (tmin:delta_t:tmax)';
        n = length(t_sig);
        image_filenames = sort(split(strip(ls([path, '/*.tif']))));   % image files (sorted)
        I = zeros(range(VerticalPixels)+1,range(HorizontalPixels)+1,n);       % matrix which array C will be converted to
        try
            for i = 1:n     %% extract pixel intensities and time delays for each frame
                image_filename_temp = erase(string(image_filenames{i}) , "'" );
                Image_handle = Tiff(image_filename_temp,'r');  % read given frame
                I(:,:,i) = read(Image_handle); % extract intensity data for given frame
            end
        catch
            status = 0;
            disp('Cannot load Tif file!')
            return;     %%% End the program if something is wrong in the file
        end
        IRF_file = strcat('IRF',num2str(param.IRF_dt),'t_',num2str(Gate),'G_',num2str(BINNING),'BIN_',param.laser,'.mat');
        if isfile(IRF_file)
            load(IRF_file,'IRF','time_axis_IRF','t_IRF');
            
            if (exist('time_axis_IRF','var'))
                t_IRF = time_axis_IRF';
            else
                t_IRF = t_IRF';
            end
        else
            disp('IRF Does not exist!')
            status = 3;
            return
        end

        if (param.BIN>1)
            new_rows = round(size(I, 1) / param.BIN);
            new_cols = round(size(I, 2) / param.BIN);
            I = squeeze(Bin(I,param.BIN,param.BIN));
            IRF = squeeze(Bin(IRF,param.BIN,param.BIN));
            HorizontalPixels = 1:new_cols;     % based on binning factor [pixels]
            VerticalPixels = 1:new_rows;         % based on binning factor [pixels]
        end
    end
    %% Performing Time Domain Analysis
    if (param.SHOW_ONLY == 0)
        I_avg = sum(I,3)/size(I,3);
        [~,idx] = max(sum(I,[1 2]));
        Peak_Img = I(:,:,idx);

        err_map = zeros(size(I,1),size(I,2)); %% Err Map is the main vector that holds the information about which pixels to fit or not
        err_map(param.Mask == 0) = 1;
        [status] = CheckTimeAxis(t_sig,t_IRF);
        if (status == 2)
            disp("Error in time axis")
            return
        end
        if (sum(err_map(:)) == 0)
            disp('Not enough pixels. Check thresholding!')
            status = 4;
            return
        end
        Obj_Mask = param.Mask;
        [status,amp_map,tau_map,err_map,chi_map,tau_seg,chi_seg,err_seg,~,I_fit_seg,I_seg_actual,tau_dist_seg,idx_vector] = TimeDomainFit(I,IRF,t_sig,t_IRF,err_map,param);
        if (status >= 5)
            disp('Errors during the fitting')
            return
        end
        time = t_sig;
        %%%% Need to check that not all tau_ map is NaN
        save(filename_saving_results,'I_avg','Peak_Img','idx_vector','param','err_seg','chi_seg','tau_map','err_map','amp_map','chi_map','HorizontalPixels','VerticalPixels','I_fit_seg','I_seg_actual','tau_seg','tau_dist_seg','time','Obj_Mask','-mat')
    end    
    %% Process the Fitted data and plotting the results 
    if(param_plot.plot_flag == 1 && param.pixelwise == 1)
        plot_data;
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
    else
        fprintf("Good time axis for signal and IRF!")
    end    
end

function binned_I = Bin(I,block_c,block_r)
    new_rows = round(size(I, 1) / block_r);
    new_cols = round(size(I, 2) / block_c);    
    binned_I = zeros(new_rows, new_cols, size(I, 3));
    for i = 1 : size(I, 3)
        binned_I(:, :, i) = blockproc(I(:, :, i), [block_r block_c], @(block_struct) mean(block_struct.data(:)));
    end
end

function [status,amp_map,tau_map,err_map,chi_map,tau_seg,chi_seg,err_seg,I_fit,I_fit_seg,I_segment_actual,tau_dist_seg,idx_vector] = TimeDomainFit(I,IRF, t_sig, t_IRF, err_map, param)
    status = 1;    %%%% No error until otherwise stated
    
    if (param.pixelwise == 1)

        amp_map = zeros(size(I,1),size(I,2),param.order);       
        chi_map = zeros(size(I,1),size(I,2));
        tau_map = zeros(size(I,1),size(I,2),param.order);    
    else
        I_fit = [];
        chi_map = [];
        tau_map = [];
        amp_map = [];
    end
    if (param.segments == 1)
        segment_map = param.Mask;
        N_segments = unique(segment_map);
        N_segments = N_segments(2:end);
        tau_dist_seg = zeros(size(N_segments,1),length(param.tau_vec));
        I_segment_actual = zeros(size(N_segments,1),length(t_sig));
        I_fit_seg = cell(size(N_segments,1),1);
        amp_seg = zeros(size(N_segments,1) ,param.order);
        chi_seg = zeros(size(N_segments,1) ,param.order);
        tau_seg = zeros(size(N_segments,1) ,param.order);
        err_seg = zeros(size(N_segments,1) ,param.order);
    else 
        tau_dist_seg = [];
        I_segment_actual = [];
        I_fit_seg  = [];
        amp_seg = [];
        chi_seg = [];
        tau_seg = [];
        err_seg = [];
    end

    if(~isempty(param.ROI))
        ROI_map = zeros(size(I,1),size(I,2));
        ROI_map(param.ROI(2):param.ROI(2)+param.ROI(4),param.ROI(1):param.ROI(1)+param.ROI(3)) = 1;
        param.PLOT_FLAG = 1;
    end
    IRF_current = squeeze(sum(IRF,[1 2]));
    
    if(param.segments == 1)     %%% Fitting segment
        for i = 1:size(N_segments,1) %% Starts from 2 since 1 is the background
            mask = double(segment_map == N_segments(i));
            % N_pixels = sum(mask(:));
            if(param.DEBUGGING == 1)
                disp('-----------------------')
                disp(strcat('Segment number: ',num2str(i),', Out of: ',num2str(size(N_segments,1))));
            end
            y_current = squeeze(sum(mask.*I,[1 2]));  %% An image
            save('Data_segment','y_current','t_sig','IRF_current','param')
            param.tau_dist = [];
            param.prior = 0;    %%% No prior model for the fitting
            
            [y_fit,amp,tau,chi,err_status] = mydeconv(y_current, t_sig, IRF_current, t_IRF,param);
            if (err_status == 0)
                I_segment_actual(i,:) = y_current;
                if (~isempty(param.tau_dist))
                    tau_dist_seg(i,:) = param.tau_dist;
                end
                I_fit_seg(i,:) = y_fit;
                chi_seg(i)  = chi;
                tau_seg(i,:) = tau;
                amp_seg(i,:)  = amp;
                err_seg(i) = param.err;
                param.err = 0;
            end
        end
    end
    if(param.pixelwise == 1)    %%% Fitting Pixles
        % tic
        if (~isempty(param.ROI))
            [rows, cols] = find(ROI_map > 0 & err_map == 0);
        else
            [rows, cols] = find(err_map == 0);
        end
        idx_vector = [rows, cols];
        N_pixels = length(rows);
        
        disp(['Number of Pixels to analyze: ' num2str(N_pixels)]);
        chi_map_vec = cell(N_pixels,1);
        amp_map_vec = cell(N_pixels,1);
        tau_map_vec = cell(N_pixels,1);
        I_fit   = cell(N_pixels,1);
        warning('off','all');
        warning('off','MATLAB:remoteparfor:ParforWorkerAborted');
        tic
        if (N_pixels == 0)
            disp('No points to analyze! Exiting!');
            status = 7;
            return;
        end
        count_bad = 0;
        for count = 1:min(N_pixels,100)
            i = rows(count);
            j = cols(count);
            y_current = squeeze(I(i,j,:));
            if (param.DEBUGGING == 1)
                mydeconv(y_current, t_sig, IRF_current, t_IRF,param);
            else
                try
                    mydeconv(y_current, t_sig, IRF_current, t_IRF,param);       
                catch
                    count_bad = count_bad + 1;
                end
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
            fprintf('ETA: from %s to %s mins. ', num2str(round(ETA/numWorkers,2)), num2str(round(ETA,2)));
        catch
            disp(strcat('ETA: ',num2str(round(ETA,2)),' mins'))
        end
        tic
        if (param.DEBUGGING == 0)
            parfor count = 1:N_pixels
                i = rows(count);
                j = cols(count);
                y_current = squeeze(I(i,j,:));      %%%% Future update: I can be stretched to a  2D vector of good pixels only to save memory
                try
                    [y_fit,amp,tau,chi,err_status] = mydeconv(y_current, t_sig, IRF_current, t_IRF, param);
                    if (err_status == 0)     %%% If we caught an error
                        I_fit{count} = y_fit;
                        amp_map_vec{count} = amp;
                        tau_map_vec{count} = tau;
                        chi_map_vec{count}  = chi;
                    end
                catch   %%% If a new unknown error happened
                    err_status = 1;
                end
                if (err_status == 1)     %%% If we caught an error
                    I_fit{count} = [];
                    amp_map_vec{count} = NaN;
                    tau_map_vec{count} = NaN;
                    chi_map_vec{count}  = NaN;
                end
            end
        else
            for count = 1:N_pixels
                i = rows(count);
                j = cols(count);
                y_current = squeeze(I(i,j,:));
                [y_fit,amp,tau,chi,err_status] = mydeconv(y_current, t_sig, IRF_current, t_IRF,param);
                if (err_status == 0)
                    I_fit{count} = y_fit;
                    amp_map_vec{count} = amp;
                    tau_map_vec{count} = tau;
                    chi_map_vec{count}  = chi;
                else     %%% If we caught an error
                    I_fit{count} = [];
                    amp_map_vec{count}  = NaN;
                    tau_map_vec{count}  = NaN;
                    chi_map_vec{count}  = NaN;
                end
            end
        end
        
        %%% Fill the output matrices
        try
            for count = 1:N_pixels
                i = rows(count);
                j = cols(count);
                chi_map(i,j)  = chi_map_vec{count};
                amp_map(i,j,:) = amp_map_vec{count};
                tau_map(i,j,:) = tau_map_vec{count};
            end
        catch
            if (param.DEBUGGING == 1)
                 keyboard;
            else
                disp('Error in copying matrix!')
                status =  7;
                return
            end
        end
        t_final = toc/60;
        disp(strcat('Done after: ',num2str(round(t_final,2)),' mins'))
    end
end

%% Extraction functions
function [freq0,delta_t,tmin,tmax,Gate,BINNING,MCP] = GetExperimentalParameters(parametersFile)
    A = fileread(parametersFile);
    A = A(~isspace(A)); % remove all spaces
    % need to extract experimental parameters
    str1 = 'TrigFreq:'; % for laser frequency
    str2 = 'Res:';  % for time steps
    str3 = 'On';    % for binning 
    str4 = 'LaserTriggerRate='; % for laser frequency (if initial identifier fails)
    str5 = 'ScanRange:'; % for laser frequency (if initial identifier fails)
    str6 = 'CombHi';
    str7 = 'MCPGain:';
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
    % extract binning specification
    indexBin_1 = strfind(A,str3);
    indexBin_2 = indexBin_1+length(str3);
    BINNING = str2double(extractBefore(A(indexBin_2:length(A)),'x'));

    indexGate = strfind(A,str6)+length(str6);
    Gate = str2double(extractBefore(A(indexGate:indexGate+5),'ps'));

    indexGate = strfind(A,str7)+length(str7);
    MCP = str2double(extractBefore(A(indexGate:indexGate+5),'V'));

end
