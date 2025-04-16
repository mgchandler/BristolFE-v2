clear;
close all;
restoredefaultpath;
addpath('../code');

%% B1_Code: Generate random grain structures with two different grain sizes

% Sample Inputs

Iteration = 1;
G1=[15]*1e-6; % Grain size 1
G2=[30]*1e-6; % Grain size 2
RUN_NUM=[1]; % Run number - see B1 for details
sample_width = 1.5e-3 /4; %
sample_depth = 6e-3 /4; % 4
sample_depth_1 = 2.5e-3 /4; %1.5
sample_depth_2 = sample_depth-sample_depth_1;
numrun = 1e6;
el_size = 2e-6; %20 wavelengths per element
TPC_r_step = 2e-6;
elementtype = 'CPE3';
sf_dr0 = [ pwd '\'];

% NOTES
% FOR 10MHz,the estimated els per wavelength are 
% delta estimate = 4.36e-5, el_size = delta/els per wavelength
% 5 els per wavelength:8-10 micrometre
% 10 els per wavelength:4-5 micrometre
% 15 els per wavelength:2.8-3.3 micrometre
% 20 els per wavelength: 2.1-2.3 micrometre
 
% Defining Material Properties
rho = 7870; %Density
col = hsv2rgb([2/3,0,0.80]); %Colour for display
name = 'Steel';
el_type = 'CPE3';
elastic_matrix = [2.0460 1.3770 2.0460 1.3770 1.3770 2.0460 0  0  0 1.2620  0       ...
    0 0  0  1.2620    0    0    0    0   0    1.2620]*1.0e+11; % stainless_steel

%SWITCHES
set_col = 1; % 1 for colored grains, 0 for grey
mat = 1; % 1 for multiple materials, 0 for one.

% Investigate the time these functions take to run. 
[mod, GRAIN_label] = fn_gen_grain_structure_and_mesh(G1,G2,RUN_NUM,elementtype,el_size,sample_width, sample_depth,sample_depth_1,sample_depth_2, mat);
[matls, dis] = fn_gen_grain_matls_2(elastic_matrix, rho, col, name, el_type, GRAIN_label, set_col, el_size);

%--------------------------------------------------------------------------
% %DEFINE THE PROBLEM

src_end_pts = [
    sample_width*0.05, sample_depth
    sample_width*0.95, sample_depth];

% %Define start of absorbing boundary region and its thickness
% abs_bdry_thickness = 7e-5;
% 
% %absorbing boundary at bottom
% abs_bdry_pts = [
%     abs_bdry_thickness, sample_depth
%     abs_bdry_thickness, abs_bdry_thickness
%     sample_width - abs_bdry_thickness, abs_bdry_thickness
%     sample_width - abs_bdry_thickness, sample_depth];
% 
% % abs_bdry_pts = [
% %     abs_bdry_thickness, sample_depth
% %     abs_bdry_thickness, 0
% %     sample_width - abs_bdry_thickness, 0 
% %     sample_width - abs_bdry_thickness, sample_depth];
% 
% mod = fn_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness);

src_dir = 2; %direction of forces applied: 1 = x, 2 = y, 3 = z (for solids)

%TEST
centre_freq = 10e6; % Centre frequency in Hz
no_cycles = 3; % Number of cycles
max_time = 4e-6; % Duration of the pulse in seconds, increase when needed
%sampling_freq = 2e9/4; % Sampling frequency in Hz

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
time_step = fn_get_suitable_time_step_2(matls, el_size);
%steps{1}.load.time = linspace(0, duration , duration * sampling_freq); %/4
steps{1}.load.time = 0: time_step/4:  max_time;
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, centre_freq, no_cycles);

% figure;
% plot(steps{1}.load.time, steps{1}.load.frcs);
% xlabel('Time (s)');
% ylabel('Amplitude');
% title('Gaussian Pulse');
% grid on;

%Also record displacement history at same points (NB there is no reason why
%these have to be same as forcing points)
steps{1}.mon.nds = steps{1}.load.frc_nds;
steps{1}.mon.dfs = steps{1}.load.frc_dfs;

%Show the mesh
figure; 
display_options.draw_elements = 1;
display_options.draw_axes = 0;
display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
display_options.node_sets_to_plot(1).col = 'r.';
h_patch = fn_show_geometry_2(mod, matls, display_options);
%saveas(gcf, 'Grains_15_30_Abs.fig')

% %--------------------------------------------------------------------------
% %RUN THE MODEL
% 
% fe_options.field_output_every_n_frames = 100;
% res = fn_BristolFE_v4(mod, matls, steps, fe_options);
% 
% % %--------------------------------------------------------------------------
% %SHOW THE RESULTS
% 
% figure;
% plot(steps{1}.load.time, sum(res{1}.dsps)); %see dsps
% xlabel('Time (s)')
% 
% plot(steps{1}.load.time, sum(res{1}.dsps));

% % Create a structure to store your data
% sig.steps_load_time = steps{1}.load.time;
% sig.sum_res_dsps = sum(res{1}.dsps);
% 
% data.mod = mod;
% data.matls = matls;
% 
% % Specify the file name to save
% % signal_filename = '01_Signal_G15_G30_F10MHz_E2.mat';
% signal_filename = sprintf('%02d_Signal_G%d_G%d_F%d_E%d.mat', Iteration, round(G1*1e6), round(G2*1e6), round(centre_freq/1e6), round(el_size*1e6));
% model_filename = sprintf('%02d_Model_Matls_G%d_G%d_F%d_E%d.mat', Iteration, round(G1*1e6), round(G2*1e6), round(centre_freq/1e6), round(el_size*1e6));
% res_filename = sprintf('%02d_Results_G%d_G%d_F%d_E%d.mat', Iteration, round(G1*1e6), round(G2*1e6), round(centre_freq/1e6), round(el_size*1e6));
% 
% folder_name = 'Grain_Boundary_Outputs';
% 
% % Specify the full paths including folder
% signal_full_path = fullfile(folder_name, signal_filename);
% model_full_path = fullfile(folder_name, model_filename);
% res_full_path = fullfile(folder_name, res_filename);
% 
% % Save the data structures to .mat files
% save(signal_full_path, 'sig');
% save(model_full_path, 'data');
% save(res_full_path, 'res');

%Animate result
% figure;
% display_options.draw_elements = 0;
% h_patch = fn_show_geometry_2(mod, matls, display_options);
% anim_options.repeat_n_times = 1;
% fn_run_animation(h_patch, res{1}.fld, anim_options);

% Mean Velocity comparison for 6mm model
% at 10e-6 element, the output velocity is:
% mean velocity at 3e-6: 
% margin of error:

