clear all;clc;format short g
% path_main = '/home/min/a/salem8/lib/Fitting/Data/2025_10_06_AmpTest_prior';
%path_main = '/local/scratch/a/salem8/Data/wenzhocopy/'; folder_name = '1256/';
% path_main = '/local/scratch/a/salem8/Data/'; folder_name = 'movethis/';
path_main = '/local/scratch/a/salem8/Data/membrane_last/';
% path_main = '/home/min/a/salem8/lib/Fitting/Data/';

% folder_name_arr = ["NoPFF_Rab5_April","PFF_Rab5_April","NoPFF_Rab5_May","PFF_Rab5_May","NoPFF_Rab7_April","PFF_Rab7_April","NoPFF_Rab7_May","PFF_Rab7_May"];
% folder_name_arr = ["Control_Rab5_April","Control_Rab5_May","Control_Rab7_April","Control_Rab5_May"];
% folder_name_arr = ["2026_03_02_021626_A53T_D5","2026_02_20_120525_A53T_D5_v2","2026_03_18_031026_Membrane_D5_PFA"];
% folder_name_arr = ["2026_02_20_120525_A53T_D5_v2"];
% folder_name_arr = ["2026_03_24_031226_A53T_D7_PFA","2026_03_12_020626_A53T_D7_PFA","2026_02_24_120525_A53T_D7_v3"];
folder_name_arr = ["2026_04_08_020626_A53T_D5_TX_manual_v2"];
% folder_name_arr = ["2026_03_24_031226_A53T_D7_PFA"];
for folder_name_idx = 1:length(folder_name_arr)
    folder_name = folder_name_arr(folder_name_idx);
    [param,param_plot] = GetIniParam();
    
    %%% For Debugging
    param.DEBUG_Count = 0;
    param.DEBUGGING = 0;        
    param.PLOT_FLAG = 0;
    
    %%% Analyzing
    param_plot.plot_AutoObjStats = 1;
    param_plot.plot_PeakLT = 1;
    param_plot.plot_PeakInt = 1;
    % param.PriorMask = [0,0,0];
    % param.FoldersToRun = [1,2,3];
    % param.Stain = 1;
    % param.Stain_Neurites = 1;
    % N_measure_arr = 88:110;
    % N_measure_arr = 17:240;
    N_measure_arr = 0;
    % N_measure_arr = [5,8];
    % param.Analysis_type  = 2;
    param.Analysis_type  = 1;
    Auto_threshold = 0;   %%% Threshold Slider if 0. Fixed Threshold at set value.
    
    %%% Thresholding
    % Run_Sim_v5(path_main,folder_name,N_measure_arr,1,Auto_threshold,param,param_plot);
    
    %% Analysis Loop
    % for N_idx = 1:length(N_measure_arr)
    %     N_measure = N_measure_arr(N_idx);
    %     Run_Sim_v5(path_main,string(folder_name),N_measure,2,0.5,param,param_plot); %% Lifetime Fitting, Analyze all
    % end

    %% Data Extraction Loop
    for N_idx = 1:length(N_measure_arr)
        N_measure = N_measure_arr(N_idx);
        Return_data_temp = Run_Sim_v5(path_main,string(folder_name),N_measure,4,0.5,param,param_plot); %% Plot Results
    end
end