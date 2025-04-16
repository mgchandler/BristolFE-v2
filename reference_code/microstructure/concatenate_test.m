el_size = 2e-6;
% Initialize variables to store data from realization 5
time_reference = [];
dsps_reference = [];

% Load data from realization 5
folder_path = 'Grain_Boundary_Outputs/E2';
filename_reference = sprintf('%02d_Signal_G15_G30_F10_E%d.mat', 5, round(el_size*1e6));
file_path_reference = fullfile(folder_path, filename_reference);

if exist(file_path_reference, 'file')
    exp_data_reference = load(file_path_reference);
    time_reference = exp_data_reference.sig.steps_load_time;
    dsps_reference = exp_data_reference.sig.sum_res_dsps;
else
    error('Data file for realization 5 does not exist.');
end

% Loop over realizations 16-20
for run_num = 16:19
    % Determine the starting index for time_reference_subset

    % Generate the file name for the current realization
    filename_current = sprintf('%02d_Signal_G15_G30_F10_E%d.mat', run_num, round(el_size*1e6));
    file_path_current = fullfile(folder_path, filename_current);
    
    % Check if the file exists
    if exist(file_path_current, 'file')
        % Load the signal data
        exp_data_current = load(file_path_current);
        time_current = exp_data_current.sig.steps_load_time;
        dsps_current = exp_data_current.sig.sum_res_dsps;
        
        % Determine the starting index for time_reference_subset
         % Determine the starting index for time_reference_subset
        start_index = numel(time_current) + 1;

        % Select the remaining portion of time_reference
        time_reference_subset = time_reference(start_index:end);
        dsps_reference_subset = dsps_reference(start_index:end);


        % Concatenate data from realization 5
        time_current = [time_current, time_reference_subset];
        dsps_current = [dsps_current, dsps_reference_subset];
        
        % Do whatever further processing you need with the concatenated data
        % For example, save it to a new file
        % Save the concatenated data
        exp_data_current.sig.steps_load_time = time_current;
        exp_data_current.sig.sum_res_dsps = dsps_current;
        % Append '_ct' to the filename before the file extension
        [~, base_filename, ext] = fileparts(filename_current);
        new_filename = [base_filename, '_ct', ext];
        new_file_path = fullfile(folder_path, new_filename);
        
        % Save the concatenated data with the new filename
        save(new_file_path, '-struct', 'exp_data_current');
    else
        warning('Data file for realization %d does not exist.', run_num);
    end
end