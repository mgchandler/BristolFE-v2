% 14/02 Testing
% Details:
% Testing stiffness matrix - work towards implementing different material
% properties as a function of grain label. 


clear;
close all;
restoredefaultpath;
addpath('../code');

%% B1_Code
%%% Generate random grain structures with two different grain sizes,
%%% and calculate the TPCF and grain size data.
%%% Then, plot TPC curve and grain size distribution.


% Defining size of grains 1 and 2
G1=[20]*1e-5;
G2=[20]*1e-5; %20

% Define run iterations 
% To be used in the incoming weeks for investigation (12/02)
RUN_NUM=[1]; 
mat = 0;

% Define bounding shape
sample_width = 3e-3;
sample_depth = 6e-3;
sample_depth_1 = 3e-3;
sample_depth_2 = 3e-3;
numrun = 1e6;

% Define element size and r_step: defined in General Notes (12/02)
el_size = 5e-5;
TPC_r_step = 3e-6;
elementtype = 'CPE3';
sf_dr0 = [ pwd '\'];

% Setting Loop from B1
for ig2=1:length(G2)
for ig1=1:length(G1)
g2 = G2(ig2); 
g1 = G1(ig1); 
for run_num=1:RUN_NUM(ig1)
    
    % Save outputs into folder for investigation
    close all; clear I0;
    filename = ['g1_' num2str(round(g1*1e6)) '_g2_' num2str(round(g2*1e6)) ...
        '_' num2str(round(sample_width*1e3)) '-' num2str(round(sample_depth_1*1e3)) '-' num2str(round(sample_depth_2*1e3))];
    sf_dr = [sf_dr0 '\OUTPUT-grainstructure\' filename '\'];
    mkdir(sf_dr);

    % correlated parameters: defines the anistropy of the grain spread. 
    grain_size_spread_1 = g1; %/5
    grain_size_spread_2 = g2; 

    % Function to build voronoi structure
    %use one for now for testing
    [VX_1,VY_1,V_1,C_1,VX_2,VY_2,V_2,C_2]=fn_grain_structure_one_grain_size...
        (g1,grain_size_spread_1,g2,grain_size_spread_2,sample_width,sample_depth,sample_depth_1,sample_depth_2);
    save([sf_dr  'grain_structure' num2str(run_num) '.mat'],...
            'VX_1','VY_1','V_1','C_1','VX_2','VY_2','V_2','C_2'); 

    % Function to generate mesh using B1. Edited to fit triangular mesh

    % Identifies if mesh is saved in output structure
    % if isfile([sf_dr0 '\OUTPUT-grainstructure\MESH\' elementtype  '_size' num2str(round(el_size*1e9)) '_' num2str(sample_width*1e6) '-' num2str(sample_depth*1e6) '.mat']);
    %     load([sf_dr0 '\OUTPUT-grainstructure\MESH\' elementtype  '_size' num2str(round(el_size*1e9)) '_' num2str(sample_width*1e6) '-' num2str(sample_depth*1e6) '.mat']);
    % else 
    % Function to generate mesh   
        if ismember(elementtype,{'CPE3', 'CPE4','CPS4','CPE4R','CPS4R'})
  
            % Function defined by Bristol_FE_v2.
            mod = fn_tri_structured_mesh(el_size,sample_width,sample_depth,mat);

            % Use naming convention from B1 to access grain structure match
            % with mesh functionality
            nodes = mod.nds;
            elements = mod.els;

            ele_cen=[(nodes(elements(:,1),1)+nodes(elements(:,2),1)+nodes(elements(:,3),1))/3,...
                (nodes(elements(:,1),2)+nodes(elements(:,2),2)+nodes(elements(:,3),2))/3];  
        end
    %     mkdir([sf_dr0 '\OUTPUT-grainstructure\MESH\']);
    %     save([sf_dr0 '\OUTPUT-grainstructure\MESH\' elementtype '_size' num2str(round(el_size*1e9)) '_' num2str(sample_width*1e6) '-' num2str(sample_depth*1e6) '.mat'],...
    %         'nodes','elements','ele_cen');
    % end

    % Assigns the mesh elements to each grain
    GRAIN_label=fn_label_grains_two_grain_size_using_VC(ele_cen,sample_depth_2,V_1,C_1,V_2,C_2);
    save([sf_dr 'grain_label' num2str(run_num) '_' elementtype '_size' num2str(round(el_size*1e9)) '.mat'],...
        'GRAIN_label');  
    
end
end
end

%--------------------------------------------------------------------------
%DEFINE THE PROBLEM

%Material properties
% Investigate suitable materials for this, as of (12/02) use steel as a
% first test

matls(1).rho = 8900; %Density
%3x3 or 6x6 stiffness matrix of material. Here it is isotropic material and
%fn_isotropic_plane_strain_stiffness_matrix(E, v) converts Young's modulus
%and Poisson's ratio into appropriate 3x3 matrix
%matls(1).D = fn_isotropic_plane_strain_stiffness_matrix(210e9, 0.3); 
elastic_matrix = [2.0460 1.3770 2.0460 1.3770 1.3770 2.0460 0  0  0 1.2620  0       ...
    0 0  0  1.2620    0    0    0    0   0    1.2620]*1.0e+11; % stainless_steel
C0=zeros(6,6);
C0([1,7,8,13,14,15,19,20,21,22,25,26,27,28,29,31,32,33,34,35,36])=elastic_matrix;
C0=C0+tril(C0',-1);
matls(1).D = C0;
matls(1).col = hsv2rgb([2/3,0,0.80]); %Colour for display
matls(1).name = 'Steel';
matls(1).el_typ = 'CPE3'; %CPE3 must be the element type for a solid

%separate into generating megastructure and ensure that i can propose the
%microstructure and use that measure grain size to make sure the structure
%is generated properly based on the ultrasound results

% use mean from grain size against PDF trend
% try different ie 20 micrometer mean and 5 micrometer standard deviation
% microstructure is that - mean grain size mean standard deviation based on
% that the single model layer will get a backscatter confrim that the
% microstrucutre is exactly what you expect - i want to create a
% microstructure propose and use code to generate and use the relationship
% to confirm that - when we get back scanning signal mathematical mode
% behind it is 
%key question - discuss how i can change the elements such that there are
%different material properties for the grains and JZ will send script
%rotation

%Define shape of model
% Bounding size defined in mesh generation for tri in B1.
% model_size = 10e-3;
% model_size = 10e-3;
% bdry_pts = [
%     0, 0 
%     model_size, 0 
%     model_size, model_size];

%Define a line along which sources will be placed to excite waves sample_width
% src_end_pts = [
%     0.3 * model_size, 0
%     0.7 * model_size, 0];

src_end_pts = [
    0, sample_depth
    sample_width, sample_depth];

src_dir = 2; %direction of forces applied: 1 = x, 2 = y, 3 = z (for solids), 4 = volumetric expansion (for fluids)

%Details of input signal
% modify based on ultrasound
% centre_freq = 5e6;
% no_cycles = 4;
% max_time = 10e-6;

centre_freq = 10e6; % Centre frequency in Hz
no_cycles = 1; % Number of cycles
duration = 5e-6; % Duration of the pulse in seconds
sampling_freq = 2e9; % Sampling frequency in Hz

%Elements per wavelength (higher = more accurate and higher computational cost)
% JZ: Use 20 elements per wavelength - W14 meeting
els_per_wavelength = 20;

%--------------------------------------------------------------------------
%PREPARE THE MESH

%Work out element size
%user defined for testing
%el_size = fn_get_suitable_el_size(matls, centre_freq, els_per_wavelength);

%Create the nodes and elements of the mesh: defined from B1
% mod = fn_tri_structured_mesh(bdry_pts, el_size);

%Identify nodes along the source line to say where the loading will be 
%when FE model is run
steps{1}.load.frc_nds = fn_find_nodes_on_line(mod.nds, src_end_pts(1, :), src_end_pts(2, :), el_size / 2);
steps{1}.load.frc_dfs = ones(size(steps{1}.load.frc_nds)) * src_dir;

%Also provide the time signal for the loading (if this is a vector, it will
%be applied at all frc_nds/frc_dfs simultaneously; alternatively it can be a matrix
%of different time signals for each frc_nds/frc_dfs
time_step = fn_get_suitable_time_step(matls, el_size);
%steps{1}.load.time = 0: time_step:  max_time;
steps{1}.load.time = linspace(0, duration, duration * sampling_freq);
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, centre_freq, no_cycles);

%Also record displacement history at same points (NB there is no reason why
%these have to be same as forcing points)
steps{1}.mon.nds = steps{1}.load.frc_nds;
steps{1}.mon.dfs = steps{1}.load.frc_dfs;

%Show the mesh
figure; 
display_options.draw_elements = 1;
display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
display_options.node_sets_to_plot(1).col = 'r.';
h_patch = fn_show_geometry(mod, matls, display_options);

%--------------------------------------------------------------------------
%RUN THE MODEL

fe_options.field_output_every_n_frames = 10;
res = fn_BristolFE_v2(mod, matls, steps, fe_options);

%--------------------------------------------------------------------------
%SHOW THE RESULTS

% Show the history output as a function of time - here we just sum over all 
% the nodes where displacments were recorded
figure;
plot(steps{1}.load.time, sum(res{1}.dsps));
xlabel('Time (s)')

%Animate result
figure;
display_options.draw_elements = 0;
h_patch = fn_show_geometry(mod, matls, display_options);
anim_options.repeat_n_times = 1;
fn_run_animation(h_patch, res{1}.fld, anim_options);

