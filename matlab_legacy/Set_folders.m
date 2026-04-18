% baseDir = '/local/scratch/a/salem8/Data/Oct25_A53T_lamp1/2025_27_10_092625_Lamp1_D2';
% baseDir = '/local/scratch/a/salem8/Data/Oct25_A53T_lamp1/2025_26_10_092625_Lamp1_D1';
% baseDir = '/local/scratch/a/salem8/Data/Oct25_A53T_lamp1/2025_25_10_041125_A53T_D4';
% baseDir = '/local/scratch/a/salem8/Data/Oct25_A53T_lamp1/2025_24_10_041125_A53T_D5_Notsure';
% baseDir = '/local/scratch/a/salem8/Data/Oct25_A53T_lamp1/2025_21_10_050225_A53T_D5_v3';
% baseDir = '/local/scratch/a/salem8/Data/Oct25_A53T_lamp1/2025_21_10_050225_A53T_D5';
% baseDir = '/local/scratch/a/salem8/Data/Oct25_A53T_lamp1/2025_27_10_092625_Lamp1_D2_v2';
% baseDir = '/local/scratch/a/salem8/Data/Control/25_20_11_070425_mt2_mvenus_bleed_thru';
% baseDir = '/local/scratch/a/salem8/Data/ForCells/Control_Rab7_April';
baseDir = '/local/scratch/a/salem8/Data/membrane_last_D7/2026_03_24_031226_A53T_D7_PFA';
% Filters you care about
filters = {'Filter2', 'Filter3', 'Filter4'};


% Get a list of all items in the base directory
items = dir(baseDir);

% Create destination folders if they don't exist
for i = 1:numel(filters)
    destDir = fullfile(baseDir, filters{i});
    if ~exist(destDir, 'dir')
        mkdir(destDir);
    end
end
% Loop through all items
for i = 1:length(items)
    name = items(i).name;
    
    % Skip '.' and '..'
    if items(i).isdir && ~strcmp(name, '.') && ~strcmp(name, '..') && ~strcmp(name,'Properties')
        % Check for 'Filter' followed by a number
        tokens = regexp(name, 'Filter(\d+)', 'tokens');
        if ~isempty(tokens)
            filterNum = tokens{1}{1};
            destFolder = fullfile(baseDir, ['Filter' filterNum]);
            
            % Only move if this filter folder exists (e.g., Filter2/3/4)
            if exist(destFolder, 'dir')
                sourcePath = fullfile(baseDir, name);
                destPath = fullfile(destFolder, name);
                
                fprintf('Moving %s → %s\n', name, destFolder);
                movefile(sourcePath, destPath);
            else
                fprintf('Skipping %s (Filter%s not in list)\n', name, filterNum);
            end
        end
    end
end

disp('Done!');
