%% Integrate all waves identical to G1, G2F, into a matrix for easy subsequent processing

% mean V for 15 elements per wavelength: 5.8416e3
% mean for 5 elements per wavelength: 5.8593e3
% mean for 20 els: 5.84e3
% v = wavelength*frequency
% wavelength = 
% Data for half size model: 
% velocity: 5.5382e3
% time peak 1: 1.99e-07
% time peak 2: 1.2756e-6
% Error 0.3%


% Define parameters
time_of_flight1 = 3e-7; % 3e-7, 1.8
time_of_flight2 = 2.1e-6; % 2.1e-6
H = 6E-3; % Example value, adjust as needed
cycles=3;
f=10e6;
lw=1.5;
el_size = 2e-6;
% Define the specific iterations you want to loop over
iter = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18];
%iter = [6, 7, 8, 12, 13, 14, 15];


% Set the total number of iterations and the number of combinations to generate
total_iters = 18;
num_combinations = 10;

% Initialize a cell array to store the combinations
combinations = cell(num_combinations, 1);

% Generate combinations
for i = 1:num_combinations
    % Generate a random permutation of the numbers from 1 to total_iters
    random_permutation = randperm(total_iters);
    % Select the first three numbers from the random permutation to form a combination
    combination = random_permutation(1:18);
    % Store the combination in the cell array
    combinations{i} = combination;
end

% Display the generated combinations
fprintf('Generated combinations:\n');
for i = 1:num_combinations
    fprintf('Combination %d: %s\n', i, mat2str(combinations{i}));
end



% Calculate the number of iterations in the specified array
iter_range = numel(iter);

% Initialize arrays with the appropriate size
V = zeros(1, iter_range);
all_data_scatter = zeros(1, iter_range);


folder_path = 'Grain_Boundary_Outputs\E2_2';
% Loop over the specified iterations

for idx = 1:iter_range
    run_num = iter(idx);
    filename = sprintf('%02d_Signal_G15_G30_F10_E%d.mat', run_num,  round(el_size*1e6));
    %filename = sprintf('%02d_SignalHalfAbs_G15_G30_F10_E%d.mat', run_num,  round(el_size*1e6)); %change when required
    file_path = fullfile(folder_path, filename);
    
    % Check if the file exists
    if exist(file_path, 'file')
        % Load the signal data
        exp_data = load(file_path);
        time = exp_data.sig.steps_load_time;
        dsps = exp_data.sig.sum_res_dsps;

        % Find the index of the maximum dsps value in the first half of the array
        max_dsps_index_1 = find(dsps == max(dsps(1:round(length(dsps)/2))));
        max_dsps_index_2 = find(dsps == max(dsps(round(length(dsps)/2):end)));
        % Find the corresponding time value
        time_peak_dsps_1 = time(max_dsps_index_1);
        time_peak_dsps_2 = time(max_dsps_index_2);
        % Calculate velocity
        V(idx) = H * 2 / (time_peak_dsps_2 - time_peak_dsps_1);

        % Find the scattering wave
        ind_scatter = find(time > time_of_flight1 & time < time_of_flight2);
        time_scatter = time(ind_scatter);
        data_scatter = dsps(ind_scatter);
        
        % Put all scattering wave in one matrix
        all_data_scatter(1:length(data_scatter), run_num) = data_scatter';
        
        %fprintf('Processed file: %s\n', filename);
    else
        fprintf('File not found: %s\n', filename);
    end
end

%% 添加滤波
igroups_V = V;
igroups_data_scatter = all_data_scatter';
%% 计算
% velocity
groups_mean_v = mean(igroups_V);
% mean
igroups_mean_data_scatter=mean(igroups_data_scatter);
% mean square
igroups_mean_square_data_scatter=mean(igroups_data_scatter.^2);
% spatial variance
igroups_spatial_variance_data_scatter=igroups_mean_square_data_scatter-igroups_mean_data_scatter.^2;
groups_spatial_variance_data_scatter=igroups_spatial_variance_data_scatter;
% mle
groups_mle_spatial_variance = fn_rlt1(igroups_spatial_variance_data_scatter');
% 计算depth_time
groups_depth_time= (time_scatter'-cycles/f/2)*groups_mean_v/2;

groups_spatial_variance_data_scatter = groups_spatial_variance_data_scatter(1:length(groups_depth_time)) ;
groups_mle_spatial_variance = groups_mle_spatial_variance(1:length(groups_depth_time))  ;

figure;
[ax, h1, h2] = plotyy(groups_depth_time*1e3,groups_spatial_variance_data_scatter,groups_depth_time*1e3,groups_mle_spatial_variance);
set(h1,'LineWidth', lw);
set(h2,'LineWidth', lw);
set(gca,'FontName','Times','FontSize',16);
xlabel('Depth (mm)','Interpreter','latex')
ylabel('Spatial variance', 'Interpreter','latex')

% Find the maximum value and its corresponding depth for the first dataset
[max_value_1, idx_1] = max(groups_spatial_variance_data_scatter);
depth_at_max_1 = groups_depth_time(idx_1) * 1e3; % Convert depth to millimeters

% Find the maximum value and its corresponding depth for the second dataset
[max_value_2, idx_2] = max(groups_mle_spatial_variance);
depth_at_max_2 = groups_depth_time(idx_2) * 1e3; % Convert depth to millimeters

% Display the results
fprintf('Peak of first dataset: %f at depth %f mm\n', max_value_1, depth_at_max_1);
fprintf('Peak of second dataset: %f at depth %f mm\n', max_value_2, depth_at_max_2);



