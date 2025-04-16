clear;
close all;
restoredefaultpath;
addpath('../code');

%% B1_Code: Generate random grain structures with two different grain sizes

% Sample Inputs

G1=[10]*1e-5; % Grain size 1
G2=[30]*1e-5; % Grain size 2
RUN_NUM=[1]; % Run number - see B1 for details
sample_width = 1e-3;
sample_depth = 2e-3;
sample_depth_1 = 1e-3;
sample_depth_2 = 1e-3;
numrun = 1e6;
el_size = 2e-5; %20 wavelengths per element
TPC_r_step = 3e-6;
elementtype = 'CPE3';
sf_dr0 = [ pwd '\'];

% Area: b = 0.00295 - 0.00291 = 4 x 10-5; h = 0.00586207 - 0.00589655 = 3.448x10-5

% Defining Material Properties
rho = 7870; %Density
col = hsv2rgb([2/3,0,0.80]); %Colour for display
name = 'Steel';
el_type = 'CPE3';
elastic_matrix = [2.0460 1.3770 2.0460 1.3770 1.3770 2.0460 0  0  0 1.2620  0       ...
    0 0  0  1.2620    0    0    0    0   0    1.2620]*1.0e+11; % stainless_steel

%SWITCHES
set_col = 0; % 1 for colored grains, 0 for grey
mat = 1; % 1 for multiple materials, 0 for one.


% Investigate the time these functions take to run. 
[mod, GRAIN_label] = fn_gen_grain_structure_and_mesh(G1,G2,RUN_NUM,elementtype,el_size,sample_width, sample_depth,sample_depth_1,sample_depth_2, mat);
[matls, dis] = fn_gen_grain_matls_2(elastic_matrix, rho, col, name, el_type, GRAIN_label, set_col, el_size);

%% TPC Testing
rows = mod.rows;
columns = mod.columns;

%-----------------------------------------------
% PDF figure test

% field_index = [dis(:).num_instances];
% grain_area = [dis(:).grain_area];
% 
% % Plot bar graph
% figure;
% bar(grain_area, field_index, 'b');
% xlabel('Grain Area');
% ylabel('Number of Grains');
% title('Grain Area vs. Field Index');
% grid on;

%--------------------------------------------------------------------------
% %DEFINE THE PROBLEM

src_end_pts = [
    sample_width*0.3, sample_depth
    sample_width*0.7, sample_depth];

%Define start of absorbing boundary region and its thickness
abs_bdry_thickness = 10e-5;
% abs_bdry_pts = [
%     abs_bdry_thickness, 0
%     sample_width - abs_bdry_thickness, 0
%     sample_width - abs_bdry_thickness, sample_depth - abs_bdry_thickness
%     abs_bdry_thickness, sample_depth - abs_bdry_thickness];

abs_bdry_pts = [
    abs_bdry_thickness, 0
    sample_width - abs_bdry_thickness, 0
    sample_width - abs_bdry_thickness, sample_depth
    abs_bdry_thickness, sample_depth];

mod = fn_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness);

src_dir = 2; %direction of forces applied: 1 = x, 2 = y, 3 = z (for solids), 4 = volumetric expansion (for fluids)

%Details of input signal
% modify based on ultrasound
% centre_freq = 6e6;
% no_cycles = 1;
% max_time = 4e-6;
%time_step = 

%TEST
centre_freq = 1e6; % Centre frequency in Hz
no_cycles = 3; % Number of cycles
duration = 4e-6; % Duration of the pulse in seconds, increase when needed
sampling_freq = 2e9/4; % Sampling frequency in Hz

% Generate time axis
%time_step = linspace(0, duration, duration * sampling_freq);

% Generate Gaussian pulse

%  The excitation was a five cycle, 190 kHz, Gaussian-windowed pulse,
% leading to an average wavelength of 27 mm for the S0 mode

%--------------------------------------------------------------------------

%Identify nodes along the source line to say where the loading will be 
%when FE model is run
steps{1}.load.frc_nds = fn_find_nodes_on_line(mod.nds, src_end_pts(1, :), src_end_pts(2, :), el_size / 2);
steps{1}.load.frc_dfs = ones(size(steps{1}.load.frc_nds)) * src_dir;

%Also provide the time signal for the loading (if this is a vector, it will
%be applied at all frc_nds/frc_dfs simultaneously; alternatively it can be a matrix
%of different time signals for each frc_nds/frc_dfs
%time_step = fn_get_suitable_time_step_2(matls, el_size);
steps{1}.load.time = linspace(0, duration , duration * sampling_freq); %/4
%steps{1}.load.time = 0: time_step:  max_time;
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, centre_freq, no_cycles);

figure;
plot(steps{1}.load.time, steps{1}.load.frcs);
xlabel('Time (s)');
ylabel('Amplitude');
title('Gaussian Pulse');
grid on;

%Also record displacement history at same points (NB there is no reason why
%these have to be same as forcing points)
steps{1}.mon.nds = steps{1}.load.frc_nds;
steps{1}.mon.dfs = steps{1}.load.frc_dfs;

%Show the mesh
figure; 
display_options.draw_elements = 1;
display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
display_options.node_sets_to_plot(1).col = 'r.';
h_patch = fn_show_geometry_2(mod, matls, display_options);

%--------------------------------------------------------------------------
%RUN THE MODEL
% change to 15 and 30 micrometre, 


fe_options.field_output_every_n_frames = 10;
res = fn_BristolFE_v4(mod, matls, steps, fe_options);

% %--------------------------------------------------------------------------
%SHOW THE RESULTS

figure;
plot(steps{1}.load.time, sum(res{1}.dsps)); %see dsps
xlabel('Time (s)')

plot(steps{1}.load.time, sum(res{1}.dsps));

% TESTING - DO NOT DELETE

% % Assuming steps{1}.load.time is your array
% time_p = steps{1}.load.time;
% res_dsps = sum(res{1}.dsps);
% percentage_to_trim = 0.02; % 5%
% 
% total_elements = numel(time_p); % Total number of elements in the array
% num_elements_to_trim = round(percentage_to_trim * total_elements); % Calculate the number of elements to trim
% trimmed_time = time_p(num_elements_to_trim + 1:end); % Trim the first num_elements_to_trim elements
% trimmed_dsps = res_dsps(num_elements_to_trim + 1:end);
% 
% processed=fn_rlt1(trimmed_dsps);
% 
% figure
% plot(trimmed_time,trimmed_dsps./max(trimmed_dsps))
% hold on
% plot(trimmed_time,processed./max(processed));
% 



% crosstalk window 
% xtalk_start=0.125e-6;
% xtalk_end=3.5e-6;
xtalk_start=0;
xtalk_end=2e-6;
% setting up a window of time within which to look for cross-talk

% load('Example_pristine_and_damage_ANTH_MLE.mat') % damanged and undamanged signal data
% load('time.mat') %time values corresponding to above
max_t=duration; % max time value to 0.0005 secs

time = steps{1}.load.time;
in_t=find(time<=max_t); 
time=time(in_t);% used to filter out time beyond a certain threshold (max time)

not_xtalk=find(time<xtalk_start | time>xtalk_end); %indices outside of crosstalk window

t_off=length(time)-length(not_xtalk)+1;

% undam=BAE_12_pristine(in_t);
% dam=BAE_12_damage_pos1(in_t);
dsps= sum(res{1}.dsps);

sub_sig=dsps(in_t); %difference between damaged and undamaged
processed=fn_rlt1(abs(sub_sig(not_xtalk)));

% figure
% plot(time,real([dam undam]))

figure
plot(time,real(sub_sig)./max(real(sub_sig)))
hold on
plot(time(t_off:end),real(sub_sig(not_xtalk))./max(real(sub_sig(not_xtalk))));
plot(time(t_off:end),processed./max(processed));




%Animate result
figure;
display_options.draw_elements = 0;
h_patch = fn_show_geometry_2(mod, matls, display_options);
anim_options.repeat_n_times = 1;
fn_run_animation(h_patch, res{1}.fld, anim_options);

