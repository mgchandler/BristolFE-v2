clearvars -except scripts_to_run
close all;
restoredefaultpath;
addpath(genpath('../code'));
% addpath(genpath('../subdoms'));

scripts_to_run = '.';

t1 = datetime("now", "Format", "hh:mm:ss.SSS");

%--------------------------------------------------------------------------
model_types = [
    "steel"
    "steel_aluminium"
    "steel_aluminium_in_hole"
    "randomly_assigned_to_grains"
    "full_anisotropic_grains"];
model_type = model_types(5);

%DEFINE THE PROBLEM
Iteration = 1;
grains = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 60, 75, 90, 120, 150, 200, 250]*1e-6;
RUN_NUM=[1]; % Run number - see B1 for details
sample_width = 1.5e-3; % 16e-3;
sample_depth = 2.5e-3; % 6.5e-3;

crack_depth = 0.5e-3;
crack_width = 400e-6;

numrun = 1e6;
% el_size = 12.5e-6; %20 elements per wavelength, 10MHz
el_size = 10e-6; %20 elements per wavelength, 5MHz
TPC_r_step = 2e-6;
elementtype = 'CPE3';
save_dir = [ pwd '\' 'models'];
if ~exist(save_dir)
    mkdir(save_dir)
end

% NOTES
% FOR 10MHz,the estimated els per wavelength are 
% delta estimate = 4.36e-5, el_size = delta/els per wavelength
% 5 els per wavelength:8-10 micrometre
% 10 els per wavelength:4-5 micrometre
% 15 els per wavelength:2.8-3.3 micrometre
% 20 els per wavelength: 2.1-2.3 micrometre
 
% Defining Material Properties
col = hsv2rgb([2/3,0,0.80]); %Colour for display

col2 = hsv2rgb([2/3,0,0.50]); %Colour for display


%SWITCHES

%--------------------------------------------------------------------------
% %DEFINE THE PROBLEM
src_dir = 2; %direction of forces applied: 1 = x, 2 = y, 3 = z (for solids)

%TEST
centre_freq = 15e6; % 15e6; % Centre frequency in Hz
fe_options.number_of_cycles = 4; % Number of cycles
max_time = 10e-6; % Duration of the pulse in seconds, increase when needed
%sampling_freq = 2e9/4; % Sampling frequency in Hz
fe_options.use_gpu_if_available = 1;

% From Pettit, good level is 1.5 * wavelength
abs_bdry_thickness = 1.5 * 6240 / centre_freq; % 2.0 * 6240 / centre_freq;

repeats = 1;
for jj = 0:repeats
    for ii = 10
        savename = [save_dir, '\titanium-backscatter-15mhz-' num2str(grains(ii)*1e6) 'umgrains-crack' num2str(crack_width*1e6) 'um-val' num2str(jj) '.mat'];
        if exist(savename, 'file')
            continue
        end
        t3 = datetime("now", "Format", "hh:mm:ss.SSS");
        clear main exp_data
        G1 = [grains(ii)]; % Grain size 1
        % Investigate the time these functions take to run. Moved from above
        switch model_type
            case "steel"
                bdry_pts = [
                    -abs_bdry_thickness, 0
                    -abs_bdry_thickness, sample_depth + abs_bdry_thickness
                    sample_width + abs_bdry_thickness, sample_depth + abs_bdry_thickness
                    sample_width + abs_bdry_thickness, 0];
                mod = fn_isometric_structured_mesh(bdry_pts, el_size);
                matls(1).rho = 8900; %Density
                matls(1).D = fn_isotropic_stiffness_matrix(210e9, 0.3);
                matls(1).name = 'Steel';
                matls(1).col = col;
                matls(1).el_typ = elementtype; %CPE3 must be the element type for a solid
                mod.el_mat_i(:) = 1;
                main.matls = matls;
            case "steel_aluminium"
                bdry_pts = [
                    -abs_bdry_thickness, 0
                    -abs_bdry_thickness, sample_depth + abs_bdry_thickness
                    sample_width + abs_bdry_thickness, sample_depth + abs_bdry_thickness
                    sample_width + abs_bdry_thickness, 0];
                mod = fn_isometric_structured_mesh(bdry_pts, el_size);
                interface_depth = 0.4 * sample_depth;
                matls(1).rho = 8900; %Density
                matls(1).D = fn_isotropic_stiffness_matrix(210e9, 0.3);
                matls(1).name = 'Steel';
                matls(1).col = col2;
                matls(1).el_typ = elementtype; %CPE3 must be the element type for a solid
                matls(2).rho = 2700; %Density
                matls(2).D = fn_isotropic_stiffness_matrix(70e9, 0.3);
                matls(2).name = 'Aluminium';
                matls(2).col = col;
                matls(2).el_typ = elementtype; %CPE3 must be the element type for a solid
                mod.el_mat_i(:) = 1;
                node_depths = mod.nds(:, 2);
                mod.el_mat_i(mean(node_depths(mod.els(:, :)), 2) < interface_depth) = 2;
                main.matls = matls;
            case "steel_aluminium_in_hole"
                bdry_pts = [
                    -abs_bdry_thickness, 0
                    -abs_bdry_thickness, sample_depth + abs_bdry_thickness
                    sample_width + abs_bdry_thickness, sample_depth + abs_bdry_thickness
                    sample_width + abs_bdry_thickness, 0];
                mod = fn_isometric_structured_mesh(bdry_pts, el_size);
                interface_depth = 0.75 * sample_depth;
                matls(1).rho = 8900; %Density
                matls(1).D = fn_isotropic_stiffness_matrix(210e9, 0.3);
                matls(1).name = 'Steel';
                matls(1).col = col2;
                matls(1).el_typ = elementtype; %CPE3 must be the element type for a solid
                matls(2).rho = 2700; %Density
                matls(2).D = fn_isotropic_stiffness_matrix(70e9, 0.3);
                matls(2).name = 'Aluminium';
                matls(2).col = col;
                matls(2).el_typ = elementtype; %CPE3 must be the element type for a solid
                mod.el_mat_i(:) = 1;
                node_depths = mod.nds(:, 2);
                mod.el_mat_i(mean(node_depths(mod.els(:, :)), 2) < interface_depth) = 2;
                main.matls = matls;
            case "full_anisotropic_grains"
                set_col = 1; % 1 for colored grains, 0 for grey
                mat = 1; % 1 for multiple materials, 0 for one.
        %         rho = 7870; %Steel Density
                rho = 4470; %Ti Density
        %         name = 'Steel';
                name = 'Titanium';
                el_type = 'CPE3';
%                 elastic_matrix = [204.60 137.70 204.60 137.70 137.70 204.60 0  0  0 126.20  0       ...
%                     0 0  0  126.20    0    0    0    0   0    126.20]*1.0e+11; % stainless_steel
                elastic_matrix = [170  92 170  70  70 192   0   0   0  52  0       ...
                    0 0  0  52    0    0    0    0   0    39]*1.0e+9; % titanium
                [mod, GRAIN_label] = fn_gen_grain_structure_and_mesh(G1,RUN_NUM,elementtype,el_size,sample_width + 2*abs_bdry_thickness, sample_depth, mat);
                mod.nds = mod.nds - [abs_bdry_thickness, 0];
                [matls, dis] = fn_gen_grain_matls_2(elastic_matrix, rho, col, name, el_type, GRAIN_label, set_col, el_size);
                %Material properties
                main.matls = cell2mat(matls);
        end
        
        
        %Define start of absorbing boundary region and its thickness
        abs_bdry_pts = [
            0, 0
            sample_width, 0
            sample_width, sample_depth
            0, sample_depth];
        
        mod = fn_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness);
        
        %Define array
        pitch = 0.21e-3;
        no_els = 2; % 64;
        array_depth = 0;
        centre = [sample_width / 2, array_depth];
        
        %Details of input signal
        time_step_safety_factor = 3;
        
        %Elements per wavelength (higher = more accurate and higher computational cost)
        
        %The default option is field_output_every_n_frames = inf, which means there
        %is no field output. Set to a finite value to get a field output. Note that
        %in subdomain models, requesting field output causes the main model to be
        %executed twice for each transducer element, once to generate the transfer
        %functions and once to generate the field output.
%         fe_options.field_output_every_n_frames = 50;
        %--------------------------------------------------------------------------
        %PREPARE THE MESH
        main.mod = mod;
        % main.matls = matls;
        
        %Timestep
        main.mod.max_safe_time_step = fn_get_suitable_time_step(main.matls(1), el_size, time_step_safety_factor);
        main.mod.design_centre_freq = centre_freq;
        time_for_full_refl = 2 * sqrt(sample_depth^2 + sample_width^2) / 6240;
        fe_options.time_pts = ceil(3.5 * (time_for_full_refl / main.mod.max_safe_time_step));
        
        %Define array
        tc = mean(1:no_els);
        for t = 1:no_els
            el_start = centre + [pitch * (t - tc - 0.5), 0];
            el_end =   centre + [pitch * (t - tc + 0.5), 0];
            [main.trans{t}.nds, s] = fn_find_nodes_on_line(main.mod.nds, el_start, el_end, el_size / 2);
            main.trans{t}.dfs = ones(size(main.trans{t}.nds)) * 2; %DF 2 is y direction
        end
        
        %Create a subdomain in the middle with a hole in surface as scatterer
        a = linspace(0, 2*pi, 361)';
        subdomain_size = 1.0e-3;
        inner_bdry = [cos(a), sin(a)] / 2 * subdomain_size + [sample_width/2, sample_depth];% .* [sample_width, sample_depth];
        scat_pts = [[-crack_width/2; crack_width/2; crack_width/2; -crack_width/2], [0; 0; -crack_depth; -crack_depth]] + [sample_width / 2, sample_depth];
        main.doms{1}.mod = fn_create_subdomain(main.mod, main.matls, inner_bdry, abs_bdry_thickness);
%         main.doms{1}.mod = fn_add_scatterer(main.doms{1}.mod, main.matls, scat_pts, 0);

        %Show the mesh
        % if ~exist('scripts_to_run') %suppress graphics when running all scripts for testing
        % 
%             figure;
%             display_options.draw_elements = 0;
%             col = 'rgbmkyc';
%             for t = 1:no_els
%                 display_options.node_sets_to_plot(t).nd = main.trans{t}.nds;
%                 display_options.node_sets_to_plot(t).col = 'r.';
%             end
%             h_patch = fn_show_geometry_2(mod, cell2mat(matls), display_options);
        % end
        %--------------------------------------------------------------------------
        
        %Run main model
        fe_options.validation_mode = false;
        main = fn_run_main_model(main, fe_options);
%         fe_options.validation_mode = true;
%         main = fn_run_main_model(main, fe_options);

%         % Reduce the main model - not all needed any longer.
%         main_reduced = main;
%         for fld =  ["nds", "els", "rows", "columns", "el_mat_i", "el_abs_i", "max_safe_time_step", "design_centre_freq"]
%             main_reduced.mod = rmfield(main_reduced.mod, fld);
%         end
%         [unique_matls, unique_matls_i] = unique(main_reduced.doms{1}.mod.el_mat_i);
%         main_reduced.matls = main_reduced.matls(unique_matls);
%         main_reduced.doms{1}.mod.el_mat_i = unique_matls_i;
%         main_reduced.mod.el_typ_i = main_reduced.mod.el_typ_i(unique_matls);
    % 
    % %Run sub-domain model
        fe_options.validation_mode = false;
        main = fn_run_subdomain_model(main, fe_options);
    
    % %Animate results if requested
    % if ~isinf(fe_options.field_output_every_n_frames)
        % figure;
        % anim_options.repeat_n_times = 1;
        % anim_options.db_range = [-40, 0];
        % anim_options.pause_value = 0.001;
        % anim_options.frame_rate = 24;
        % % anim_options.mp4_out = './microstructure_6.5mm_8el_5mhz.mp4';
        % h_patches = fn_show_geometry_with_subdomains(main, anim_options);
        % fn_run_subdomain_animations(main, h_patches, anim_options);
    % end
    % 
    % %Run validation model
        fe_options.validation_mode = 1;
        main = fn_run_main_model(main, fe_options);
    
    %     save([save_dir, '\titanium_halfspace_5mhz_pristine_' num2str(ii) '.mat'], "main")
    % 
    % %Animate validation results if requested
    % if ~exist('scripts_to_run') %suppress graphics when running all scripts for testing
        % if ~isinf(fe_options.field_output_every_n_frames)
            % figure;
        %     anim_options.repeat_n_times = 1;
        %     anim_options.db_range = [-60, 0];
        %     anim_options.pause_value = 0.001;
        % anim_options.mp4_out = './microstructure_validation_6.5mm_8el_5mhz.mp4';
        %     h_patches = fn_show_geometry(main.doms{1}.val_mod, main.matls, anim_options);
        %     for t = 1:no_els
        %         fn_run_animation(h_patches, main.doms{1}.val.trans{t}.fld, anim_options);
        %     end
        % end
    % end
    % 
    % %View the time domain data and compare wih validation
        % figure;
        % i = max(find(abs(main.inp.sig) > max(abs(main.inp.sig)) / 1000));
        % mv = max(abs(sum(main.doms{1}.res.fmc.time_data(i:end,: ), 2)));
        % plot(main.doms{1}.res.fmc.time, sum(main.doms{1}.res.fmc.time_data, 2) / mv, 'k');
        % hold on;
        % plot(main.doms{1}.val.fmc.time, sum(main.doms{1}.val.fmc.time_data, 2) / mv, 'b');
        % plot(main.doms{1}.res.fmc.time, (sum(main.doms{1}.res.fmc.time_data, 2) - sum(main.doms{1}.val.fmc.time_data, 2)) / mv, 'r');
        % ylim([-1,1]);
        % legend('Sub-domain method', 'Validation', 'Difference');
        % savefig([save_dir '/Backscattering-abslayer-' num2str(abs_bdry_thickness / 6240 * centre_freq) 'lambda-' num2str(centre_freq*1e-6) 'mhz-titanium.fig'])
    
        exp_data = main.doms{1}.val.fmc;
        exp_data.ph_velocity = 6240;
%         save(savename, 'exp_data')
        t4 = datetime("now", "Format", "hh:mm:ss.SSS");
        disp(t4 - t3)
    end
end

t2 = datetime("now", "Format", "hh:mm:ss.SSS");

disp(t2 - t1)