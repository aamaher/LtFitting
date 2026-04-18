clear all;clc;format short g
% path_main = '/home/min/a/salem8/lib/Fitting/Data/2025_10_06_AmpTest_prior';
path_main = '/local/scratch/a/salem8/Data/membrane_last'; 
% folder_name_arr = ["NoPFF_Rab5_April","NoPFF_Rab5_May","NoPFF_Rab5_May","NoPFF_Rab7_April","PFF_Rab5_April","PFF_Rab5_May","PFF_Rab7_April"];
% folder_name_arr = ["NoPFF_Rab5_April","PFF_Rab5_April","NoPFF_Rab5_May","PFF_Rab5_May","NoPFF_Rab7_April","PFF_Rab7_April","NoPFF_Rab7_May","PFF_Rab7_May"];
folder_name_arr = ["2026_02_20_120525_A53T_D5_v2","2026_03_02_021626_A53T_D5","2026_03_18_031026_Membrane_D5_PFA"];
% folder_name_arr = ["NoPFF_Rab5_April","PFF_Rab5_April","Control_Rab5_unknown"];
% folder_name_arr = ["NoPFF_Rab5_May","PFF_Rab5_May"];
% folder_name_arr = ["Control_Rab5_unknown"];
% folder_name_arr = ["NoPFF_Rab7_May"];

% folder_name_arr = ["NoPFF_Rab7_May"];
for folder_name_idx = 1:length(folder_name_arr)
    folder_name = folder_name_arr(folder_name_idx);
    [param,param_plot] = GetIniParam();  
    param.Stain = 1;
    param.Stain_Neurites = 1;
    % param.Rerun_Mask = 1;
    param.Rerun = 1;
    % param_plot.plot_Colocalize = 1;
    %%% For Debugging
    param.DEBUG_Count = 0;
    param.DEBUGGING = 0;        
    param.PLOT_FLAG = 0;
    %%% Analyzing
    param_plot.plot_AutoObjStats = 1;
    param.PriorMask = [0,0,0];
    param.FoldersToRun = [1,2,3];
    N_measure_arr = 1:120;
    param.Analysis_type  = 3;
    Auto_threshold = 0;   %%% Threshold Slider if 0. Fixed Threshold at set value.
    
    %%% Thresholding
    % Run_Sim_v5(path_main,folder_name,N_measure_arr,1,Auto_threshold,param,param_plot);
    
    Return_pixel_flag = 0;
    Return_EM_flag = 0;
    for N_idx = 1:length(N_measure_arr)
        N_measure = N_measure_arr(N_idx);
        %%%%% Pixel-Wise
        % param.EM_Algorithm = 0;
        % % Run_Sim_v5(path_main,string(folder_name),N_measure,2,0.5,param,param_plot); %% Lifetime Fitting, Analyze all
        % Return_data_temp = Run_Sim_v5(path_main,string(folder_name),N_measure,4,0.5,param,param_plot); %% Plot Results
        % Return_data_temp = Return_data_temp{1};
        % if (isfield(Return_data_temp,'theta_avg'))
        %     if (Return_pixel_flag == 0)
        %         fields = fieldnames(Return_data_temp);
        %         Return_data_pixel{folder_name_idx} = Return_data_temp;
        %         Return_pixel_flag = 1;
        %     elseif (Return_pixel_flag == 1)
        %         if numel(fields) < numel(fieldnames(Return_data_temp))
        %             fields = unique(cat(1,fields,fieldnames(Return_data_temp)));
        %         end
        %         for idx_field = 1:size(fields,1)
        %             field  = fields(idx_field);
        %             if strcmp(field,'status')
        %                 continue;
        %             end
        %             field_value = getfield(Return_data_pixel{folder_name_idx},field{1});
        %             added_value = getfield(Return_data_temp,field{1});
        %             new_field_value = cat(1,field_value,added_value);
        %             Return_data_pixel{folder_name_idx} = setfield(Return_data_pixel{folder_name_idx},field{1},new_field_value);
        %         end
        %     end
        % end
        %%%%% EM-Algorithm
        param.EM_Algorithm = 1;
        Run_Sim_v5(path_main,string(folder_name),N_measure,2,0.5,param,param_plot); %% Lifetime Fitting, Analyze all
        Return_data_temp = Run_Sim_v5(path_main,string(folder_name),N_measure,4,0.5,param,param_plot); %% Plot Results
        Return_data_temp = Return_data_temp{1};
        if (isfield(Return_data_temp,'theta_avg'))
            if (Return_EM_flag == 0)
                fields = fieldnames(Return_data_temp);
                Return_data_EM{folder_name_idx} = Return_data_temp;
                Return_EM_flag = 1;
            elseif (Return_pixel_flag == 1)
                for idx_field = 1:size(fields,1)
                   field  = fields(idx_field);
                   if strcmp(field,'status')
                       continue;
                   end
                   field_value = getfield(Return_data_EM{folder_name_idx},field{1});
                   added_value = getfield(Return_data_temp,field{1});
                   new_field_value = cat(1,field_value,added_value);
                   Return_data_EM{folder_name_idx} = setfield(Return_data_EM{folder_name_idx},field{1},new_field_value);
                end
            end    
        end
    end
    % Return_data_pixel{folder_name_idx}.Type = "Exp";
    % Return_data_pixel{folder_name_idx}.str_leg = strcat("",folder_name_arr(folder_name_idx));
    Return_data_EM{folder_name_idx}.Type = "Exp";
    Return_data_EM{folder_name_idx}.str_leg = strcat("",folder_name_arr(folder_name_idx));
end

%%
save("Exp_data_EM.mat","Return_data_EM","folder_name_arr","-mat")
plotPair (Return_data_EM{1},Return_data_EM{2},"4_Rab5");
plotPair (Return_data_EM{3},Return_data_EM{4},"5_Rab5");
plotPair (Return_data_EM{5},Return_data_EM{6},"4_Rab7");
plotPair (Return_data_EM{7},Return_data_EM{8},"5_Rab7");
function plotPair(a1,a2,str_name)
    % plot_FRET_boxes_cells_twoPanel({a1,a2})
    % fig = gcf;
    % fig.Position(3) = 1362;
    % fig.Position(4) = 350;
    % savefig(fig, 'EM_Mem_Experimental_new.fig');
    plot_FRET_boxes_cells_v2({a1,a2})
    fig = gcf;
    fig.Position(3) = 1362;
    fig.Position(4) = 550;
    savefig(fig, strcat("EM_Mem_Experimental_",str_name,".fig"));
end