time_of_flight1 = 4e-7; % 3e-7
time_of_flight2 = 2.1e-6; % 2.1e-6
H = 6E-3; % Example value, adjust as needed
cycles=3;
f=10e6;
lw=1.5;
el_size = 2e-6;
run_num = 6;

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
plot(time*1e6, dsps*1e6, 'LineWidth', 1.3);
hold on;
plot(time*1e6, hilbert_env*1e6, 'LineWidth', 1.3);
xline(0.3, '--k', 'LineWidth', 1.5); % Add vertical dashed line at x = 0.4, with black color
xline(2.1, '--k', 'LineWidth', 1.5); % Add vertical dashed line at x = 2.1, with black color
text(0.8, 0.05, 'Backscatter Signal' ,'Interpreter', 'latex','FontSize',14);
text(0.4, -0.07, 'Front Wall' ,'Interpreter', 'latex','FontSize',14);
text(1.5, -0.07, 'Back Wall' ,'Interpreter', 'latex','FontSize',14);

xlabel('Time');
ylabel('Amplitude');

set(gca,'FontName','Times','FontSize',14);

xlabel('Time ($\mu$s)','Interpreter', 'latex');
ylabel('Amplitude ($\mu$m)','Interpreter', 'latex');
xlim([0,2.6])
legend('Received Signal', 'Hilbert Envelope', 'Interpreter', 'latex','FontSize',14)

fig3 = gcf; % Get current figure handle
% current_position = fig3.Position; % Get current figure position
% fig3.Position = [current_position(1:2), current_position(3)*2, current_position(4)]; % Adjust figure position

    % Save the figure
saveas(fig3, 'Signal_Short.png');

% % Save the data into a structure
% data2.time = time;
% data2.dsps = dsps;
% data2.hilbert_env = hilbert_env;
% 
% % Specify the filename
% filename = 'data2.mat';
% 
% % Save the data
% save(filename, 'data2');
