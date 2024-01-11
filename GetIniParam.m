function [param,param_plot] = GetIniParam(~,~)
    param = struct( ...
... %%% Simulation flags
        'DEBUGGING',0, ...              %%% Flag used for printing fitting values and plotting
        'PLOT_FLAG',0, ...              %%% Flag used to allow plotting at each fit
        'PLOT_Freq',0, ...              %%% Flag used During Debug, plot the frequency domain signal
        'PLOT_Tau',0, ...               %%% Flag used During Debug, plot the frequency domain signal
        'testtau',0, ...                %%% Unavailable now, runs fitting 1000 times and draw histogram
        'SHOW_ONLY',0, ...              %%% For internal use only. Flag used for skipping analysis and saving of fig
... %%% Parameters for processing the data
        'THRESHOLD',[], ...             %%% Threshold value used. Only used in saving file name since the mask is used for filtering
        'BIN',1, ...                    %%% Used if extra binning is allowed
        'ROI',[], ...                   %%% Region of interest for fit.
        'Mask',[], ...                  %%% Mask used for fitting. Will ignore fitting anything that is not in the mask. Used to identify objects
... %%% main Fitting parameters
        'segments',0, ...               %%% Flag used to do segment fitting, used as a prior if pixelwise is also active
        'pixelwise',1, ...              %%% Flag used to do pixel wise fitting
        'sigma_max',0.05, ...           %%% Max allowed standard deviation
        'mainlb',0.2, ...               %%% Lifetime lower bound
        'mainub',7, ...                 %%% Lifetime ubber bound
        'mainIC',2, ...                 %%% Lifetime Initial condition
        'order',1, ...                  %%% Fitting order
        'tau_vec',0.1:0.01:10, ...      %%% Lifetime vector used for plotting distribution, affect speed significantly in continouous fitting
... %%% Parameters to be used during simulation, can be separated for efficieny
        'lb',0.2, ...                   %%% Lifetime lower bound temp for when a prior is used
        'ub',10, ...                    %%% Lifetime ubber bound temp for when a prior is used
        'IC',2.5, ...                   %%% Lifetime Initial condition for when a prior is used
        'err',0, ...                    %%% Used to mark if there is an error in the current fit        
        't_shift_max', 0.25,...         %%% Maximum range for shif tin bounds
        'amp_ratio_max', 0.15,...       %%% Maximum range for shif tin bounds
        'y_shift_max', 5,...         %%% Maximum range for shif tin bounds
        'tau_dist',[], ...              %%% The current lifetime distribution. Used in both discrete and continouous
        'Trials',1, ...                 %%% Not used, number of times the fitting is ran until good fit is reached
        'prior',0, ...                  %%% Decide whether a prior is used or not
... %%%%% ...
        'global',0,...                  %%% Global or local analysis during fitting
        'parallel',1, ...               %%% Use parllel processing, currently only 1 available
        'Method',2, ...                 %%% Fitting method: 1) LMS, 2) Gaussian MLE with thermal noise, 3) MLE Poission, 4) MLE Gauss
        'cont',0, ...                   %%% Continuous or discrete lifetime fitting
        'cont_type',1, ...              %%% Type of continious lifetime distribution used 1) Gaussian, 2) Gamma, 3) Inverse gamma, 4) Log-normal
        'laser','SP',...                %%% Laser used for selection of appropriate IRF, write 'SP' or 'NKT'
        'mcp',[],...                    %%% Automatically read from the file and replaced during the running of the code
        'IRF_dt',50);                   %%% Set the dt of the IRF you want to use
    
    %%% Final plot parameters
    param_plot = struct( ...
        'plot_flag',1,...               %%% Flag whether to report or not
        'figformat','.png',...          %%% 1) .fig, 2) .png
        'SaveInFolder',0,...            %%% For saving results plots in the respective folder name similar to the data 1) yes, 0) no
        'hbins_tau', 0.1:0.01:10, ...   %%% Histogram bins used for calculating the results
        'plot_mix',0, ...               %%% Plot all All info
        'plot_hist',0, ...              %%% Plot Histogram
        'plot_chi',1, ...               %%% Plot Chi
        'plot_PeakInt',0, ...           %%% Plot Peak Intensity image
        'plot_PeakLT',0, ...            %%% Plot Peak intensity and Lifetime overlay
        'Stack_figure',0,...            %%% Only in plotting Histogram Images
        'scale_bar',20,...              %%% length of scale bar in micrometers
        'scale_bar_Linewidth',5,...     %%% Scale bar linewidth
        'Magnification',40, ...         %%% Magnification used, for plotting scale bar.
        'Fit_Histogram',0,...           %%% Fit the histogram to gaussian on the plot
        'Crop_plot',0,...               %%% Cropping the output plot, 0) No crop. 1) Manual Crop (Not set yet) 2) Automatic crop
        'fontsz',40);                   %%% Font size used in plots
end

%%% DEBUGGING flag: enables debugging mode if DEBUGGING = 1
%%% In debugging mode: 
% 1) No try-catch statement, errors will show
% 2) Fitted lifetime will be plotted for each point and the simulation will
% pause afterwards
% 3) Only the Region of interest will be analyzed. You can use this mode to
% pinpoint a specific group of pixels for debugging
% ROI WORKS AS FOLLOWS: [ymin, ymax, xmin, xmax]
% keyboard