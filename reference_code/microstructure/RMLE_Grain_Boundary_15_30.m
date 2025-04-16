close all
clear all


% Load results
load('Grain_Boundary_15_30_10MHz_20Els.mat') % Load damaged and undamaged signal data




full_time = data.steps_load_time;
signal = data.sum_res_dsps;

min_t = 0.3e-6; % Remove peak caused by input pulse
max_t = 2.5e-6; % Estimated time after first backscatter

% Find indices of time points within the specified range
in_t = find(full_time > min_t & full_time <= max_t); 
time = full_time(in_t);

processed=fn_rlt1(abs(signal(in_t)));

figure
plot(full_time,real(signal)./max(real(signal)))
hold on
plot(time,real(signal(in_t))./max(real(signal(in_t))));
plot(time,processed./max(processed));

