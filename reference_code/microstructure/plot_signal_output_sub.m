time_of_flight1 = 4e-7; % 3e-7
time_of_flight2 = 2.1e-6; % 2.1e-6
H = 6E-3; % Example value, adjust as needed
cycles=3;
f=10e6;
lw=1.5;
el_size = 2e-6;
run_num = 7;

% Define the path to the folder containing the signal data
folder_path = 'Grain_Boundary_Outputs/E2';

% Loop over 20 runs

    % Generate the file name based on the naming convention
%filename = sprintf('%02d_Signal_G15_G30_F10_E%d_ct.mat', run_num,  round(el_size*1e6));
filename = sprintf('%02d_Signal_G15_G30_F10_E%d.mat', run_num,  round(el_size*1e6)); %change when required
file_path = fullfile(folder_path, filename);
    
    % Check if the file exists
    if exist(file_path, 'file')
        % Load the signal data
        exp_data = load(file_path);
        time = exp_data.sig.steps_load_time;
        dsps = exp_data.sig.sum_res_dsps;

        sub = 0.025 * dsps;
        dsps = dsps - sub;

        max_dsps_index_1 = find(dsps == max(dsps(1:round(length(dsps)/2))));
        max_dsps_index_2 = find(dsps == max(dsps(round(length(dsps)/2):end)));

        max_dsps_wave_1 = find(dsps == max(dsps(1:round(length(dsps)/32))));
        max_dsps_wave_2 = find(dsps == max(dsps(round(length(dsps)/32):end)));
        % Find the corresponding time value
        time_peak_dsps_1 = time(max_dsps_index_1);
        time_peak_dsps_2 = time(max_dsps_index_2);

        time_peak_wave_1 = time(max_dsps_wave_1);
        time_peak_wave_2 = time(max_dsps_wave_2);
        % Calculate velocity
        V = H * 2 / (time_peak_dsps_2 - time_peak_dsps_1);
        wavelength = (time_peak_wave_2 - time_peak_wave_1) * V ;
        % Calculate the Hilbert envelope
        analytic_signal = hilbert(dsps);
        hilbert_env = abs(analytic_signal);
        

        max_dsps_hilbert_index_1 = find(hilbert_env == max(hilbert_env(1:round(length(hilbert_env)/2))));
        max_dsps_hilbert_1 = hilbert_env(max_dsps_hilbert_index_1);
        time_hilbert_1 = time(max_dsps_hilbert_index_1);

        max_dsps_hilbert_index_2 = find(hilbert_env == max(hilbert_env(round(length(hilbert_env)/2):end)));
        max_dsps_hilbert_2 = hilbert_env(max_dsps_hilbert_index_2);
        time_hilbert_2 = time(max_dsps_hilbert_index_2);

        hilbert_v = H * 2 / (time_hilbert_2 - time_hilbert_1);
        hilbert_diff = max_dsps_hilbert_1 - max_dsps_hilbert_2;
        
    end


figure;
plot(time, dsps, 'b', 'LineWidth', 1);
hold on;
plot(time, hilbert_env, 'r', 'LineWidth', 1);
xlabel('Time');
ylabel('Amplitude');



% Save the data into a structure
data23.time = time;
data23.dsps = dsps;
data23.hilbert_env = hilbert_env;

% Specify the filename
filename = 'data23.mat';

% Save the data
save(filename, 'data23');
