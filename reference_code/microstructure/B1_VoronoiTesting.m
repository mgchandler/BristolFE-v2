%%% Generate random grain structures with two different grain sizes,
%%% and calculate the TPCF and grain size data.
%%% Then, plot TPC curve and grain size distribution.
%%%%%%%%%
% close all;
clear;
G1=[15]*1e-6;
G2=[30]*1e-6; %leave at a larger value
%RUN_NUM=[100,100,100,100,100,100];
RUN_NUM = [1];
sample_width = 1e-4; %
sample_depth = 4e-4; %optimum 3.5
sample_depth_1 = 1.5e-4; %optimum 1.3
sample_depth_2 = sample_depth-sample_depth_1;
numrun = 1e6;
element_size = 2e-6;
TPC_r_step = 1e-6;
elementtype = 'CPE4';
sf_dr0 = [ pwd '\'];

% for sample depth, do equal number of grains per element
%%
for ig2=1:length(G2)
for ig1=1:length(G1)
g2 = G2(ig2); 
g1 = G1(ig1); 
for run_num=1:RUN_NUM(ig1)
    close all; clear I0;
    filename = ['g1_' num2str(round(g1*1e6)) '_g2_' num2str(round(g2*1e6)) ...
        '_' num2str(round(sample_width*1e3)) '-' num2str(round(sample_depth_1*1e3)) '-' num2str(round(sample_depth_2*1e3))];
    sf_dr = [sf_dr0 '\OUTPUT-grainstructure\' filename '\'];
    mkdir(sf_dr); 
    %% correlated parameters 
    grain_size_spread_1 = g1;
    grain_size_spread_2 = g2; 

    % Why is this divided by 5?
    %% create new set of random grains
    [VX_1,VY_1,V_1,C_1,VX_2,VY_2,V_2,C_2]=fn_grain_structure_two_grain_size_3...
        (g1,grain_size_spread_1,g2,grain_size_spread_2,sample_width,sample_depth,sample_depth_1,sample_depth_2);
    save([sf_dr  'grain_structure' num2str(run_num) '.mat'],...
            'VX_1','VY_1','V_1','C_1','VX_2','VY_2','V_2','C_2'); % square elements
%     %% generate mesh
%     % if isfile([sf_dr0 '\OUTPUT-grainstructure\MESH\' elementtype  '_size' num2str(round(element_size*1e9)) '_' num2str(sample_width*1e6) '-' num2str(sample_depth*1e6) '.mat']);
%     %     load([sf_dr0 '\OUTPUT-grainstructure\MESH\' elementtype  '_size' num2str(round(element_size*1e9)) '_' num2str(sample_width*1e6) '-' num2str(sample_depth*1e6) '.mat']);
%     % else
%         if ismember(elementtype,{'CPE4','CPS4','CPE4R','CPS4R'}) %element types from abaqus (defined by abaqus)
%             [nodes,elements,columns,rows] = fn_make_quad_mesh(element_size,sample_width,sample_depth);
%             ele_cen=[(nodes(elements(:,1),1)+nodes(elements(:,2),1)+nodes(elements(:,3),1)+nodes(elements(:,4),1))/4,...
%                 (nodes(elements(:,1),2)+nodes(elements(:,2),2)+nodes(elements(:,3),2)+nodes(elements(:,4),2))/4]; % if triangle on brisol, change fn_make_quad into make triangle. 
%         end
%         mkdir([sf_dr0 '\OUTPUT-grainstructure\MESH\']);
%         save([sf_dr0 '\OUTPUT-grainstructure\MESH\' elementtype '_size' num2str(round(element_size*1e9)) '_' num2str(sample_width*1e6) '-' num2str(sample_depth*1e6) '.mat'],...
%             'nodes','elements','ele_cen','columns','rows'); %use this mesh in bristolFE, either use grain and use bristol fe to mesh or use this mesh directly
%     %end
%     %% Assigns the mesh elements to each grain
%     GRAIN_label=fn_label_grains_two_grain_size_using_VC(ele_cen,sample_depth_2,V_1,C_1,V_2,C_2);
%     save([sf_dr 'grain_label' num2str(run_num) '_' elementtype '_size' num2str(round(element_size*1e9)) '.mat'],...
%         'GRAIN_label');   % so theres elements within each grain (grain can be large so there are numerous lemeents in the mesh), on bristolfe, assign properties to run the model.
%     % for bristolFE, what type of element does it support? either quadratic
%     % of square generated here 
%     %% Material properties to grain label
%     elastic_matrix = [2.0460 1.3770 2.0460 1.3770 1.3770 2.0460 0  0  0 1.2620  0       ...
%     0 0  0  1.2620    0    0    0    0   0    1.2620]*1.0e+11; % stainless_steel
%     density = 7870;
% 
%     C0=zeros(6,6);
%     C0([1,7,8,13,14,15,19,20,21,22,25,26,27,28,29,31,32,33,34,35,36])=elastic_matrix;
%     C0=C0+tril(C0',-1);
% 
%     %%%% euler methods, orientation angles have an uniform distribution
%     % rand('seed',sum(230*clock));
%     rand('seed',1e6);
%     euler_angle(:,1)=2*pi*rand(length(GRAIN_label),1);
%     euler_angle(:,2)=acos(2*rand(length(GRAIN_label),1)-1);
%     euler_angle(:,3)=2*pi*rand(length(GRAIN_label),1);
% 
%     CA=fn_euler_c(C0,euler_angle(:,1),euler_angle(:,2),euler_angle(:,3));
% 
%     for i=1:length(GRAIN_label)
%         % temp_material.material{i}.paramsType=anisotropic;
%         temp_material.material{i}.paramValues=[CA{i}([1,7,8,13,14,15,19,20,21,22,25,26,27,28,29,31,32,33,34,35,36]) density];
%         temp_material.matTypeRefs(GRAIN_label{i})=i;
%     end
%     %% convert image data 
%     rows_1 = round(rows*sample_depth_1/sample_depth_2);
%     for i=1:length(GRAIN_label)
%         I0(GRAIN_label{i})=i;
%     end
%     I1=reshape(I0,columns,rows);
%     pixels_size = sample_depth/rows;%% 1 pixels 
%     %% Measure TPC points for GRAIN SIZE 1
%     I2_1 = I1(:,rows_1+1:end);
%     R_1 = [TPC_r_step:TPC_r_step:g1*2]; 
%     r_1 = unique(round(R_1/pixels_size)); % unit: pixels   
%     [r_1, wr_1] = fn_calcTPCImage(I2_1, r_1, numrun);
%     r_1  = [0,r_1];
%     wr_1 = [1,wr_1];
%     %% Measure TPC points for GRAIN SIZE 2
%     I2_2 = I1(:,1:rows_1);
%     R_2 = [TPC_r_step:TPC_r_step:g2*2];
%     r_2 = unique(round(R_2/pixels_size)); % unit: pixels
%     [r_2, wr_2] = fn_calcTPCImage(I2_2, r_2, numrun);
%     r_2  = [0,r_2];
%     wr_2 = [1,wr_2];
%     %% Calculate the correlation length from the TPC curve
%     ind_1 = find(wr_1 >= 0.1);
%     aL_1 = interp1(wr_1(ind_1), r_1(ind_1)*pixels_size, exp(-1));
%     ind_2 = find(wr_2 >= 0.1);
%     aL_2 = interp1(wr_2(ind_2), r_2(ind_2)*pixels_size, exp(-1));
%     %% Measure grain size distribution for GRAIN SIZE 1
%     [grain_size_1] = fn_calcGrainSizeImage(I2_1);
%     grain_size_1 = grain_size_1*sqrt((sample_depth/rows)*(sample_width/columns));
%     %% Measure grain size distribution for GRAIN SIZE 2
%     [grain_size_2] = fn_calcGrainSizeImage(I2_2);
%     grain_size_2 = grain_size_2*sqrt((sample_depth/rows)*(sample_width/columns));
%     %% Calculate the mean grainsize 
%     ad_1 = mean(grain_size_1*1e6);
%     ad_2 = mean(grain_size_2*1e6);
% 
%     std_1 = std(grain_size_1*1e6);
%     std_2 = std(grain_size_2*1e6);
% 
%     %% PlotTPC curve
%     % fwave = figure('Visible','on');
%     % plot(r_1*pixels_size*1e6, wr_1,'Linewidth',2); hold on;
%     % plot(r_2*pixels_size*1e6, wr_2,'Linewidth',2);
%     % xlabel('Separation Distance (\mum)');
%     % ylabel('Spatial Correlation Function');
%     % box on;
%     % set(gca,'FontName','Times','FontSize',16);
%     % title(['g1-' num2str(g1*1E6) ' g2-' num2str(g2*1E6) ' ' num2str(round(sample_width*1e3))...
%     %     '-' num2str(round(sample_depth_1*1e3)) '-' num2str(round(sample_depth_2*1e3)) ' ' num2str(run_num)]);
%     % legend(['g1-' num2str(g1*1E6) '  \muL1-' num2str(round(aL_1*1E6,1))],...
%     %     ['g2-' num2str(g2*1E6) '  \muL2-' num2str(round(aL_2*1E6,1))],'location','best');
%     % print(fwave, [sf_dr '\TPC' num2str(run_num) '.emf'], '-dmeta');
% 
%     % Create a figure with LaTeX font
%     fwave = figure('Visible','on');
%     set(fwave, 'DefaultTextInterpreter', 'latex');
% 
%     % Plot the data
%     plot(r_1*pixels_size*1e6, wr_1,'Linewidth',2, 'Color', [0, 0.25, 0.5]); hold on;
%     plot(r_2*pixels_size*1e6, wr_2,'Linewidth',2, 'Color',[0.6, 0.6, 0.6]);
% 
%     % Set labels with LaTeX font
%     xlabel('Separation Distance ($\mu$m)', 'Interpreter', 'latex');
%     ylabel('Spatial Correlation Function', 'Interpreter', 'latex');
% 
%     % Customize figure
%     box on;
%     set(gca,'FontName','Times','FontSize',16);
%     title(['$g_1-' num2str(g1*1E6) '$, $g_2-' num2str(g2*1E6) '$, '...
%         num2str(round(sample_width*1e3)) '-' num2str(round(sample_depth_1*1e3)) '-'...
%         num2str(round(sample_depth_2*1e3)) ' ' num2str(run_num)], 'Interpreter', 'latex');
% 
%     % Add legend with LaTeX font
%     legend(['GS1: $ \ \mu_1 =' num2str(round(ad_1,1)) '\mathrm{\mu m},  L_1 =' num2str(round(aL_1*1E6,1)) '\mathrm{\mu m}$'],...
%         ['GS2: $ \ \mu_2 =' num2str(round(ad_2,1)) '\mathrm{\mu m},  L_2 =' num2str(round(aL_2*1E6,1)) '\mathrm{\mu m}$'],'location','northeast', 'Interpreter', 'latex',  'FontSize', 12);
% 
%     % Save figure with LaTeX font
%     print(fwave, [sf_dr '\TPC' num2str(run_num) '.emf'], '-dmeta');
% 
%     %% Plot grain size distribution
% 
%     % Generate x-values for plotting the bell curve
%     x_min = min([min(grain_size_1), min(grain_size_2)]) * 1e6;
%     x_max = max([max(grain_size_1), max(grain_size_2)]) * 1e6;
%     x_values = linspace(x_min, x_max, 1000);
% 
%     % Compute the bell curve using the normal distribution formula
%     bell_curve1 = normpdf(x_values, ad_1, std_1);
%     bell_curve2 = normpdf(x_values, ad_2, std_2);
% 
%     fwave = figure('Visible','on');
%     histogram(grain_size_1*1e6,'binwidth',1,'Normalization', 'pdf', 'FaceColor', [0, 0.25, 0.5]); hold on; %'FaceColor', [0, 0.2, 0.5]
%     hold on
%     histogram(grain_size_2*1e6,'binwidth',1,'Normalization', 'pdf', 'FaceColor',[0.7, 0.7, 0.7]); %'FaceColor',[0.5, 0.5, 0.5]
%     plot(x_values, bell_curve1, 'LineWidth', 1, 'Color', 'r' );
%     plot(x_values, bell_curve2, 'LineWidth', 1, 'Color', 'r' );
% 
%     xlabel('Grain Size ($\mathrm{\mu m}$)','Interpreter','latex');
%     ylabel('Probability Density Function','Interpreter','latex');
%     box on;
%     set(gca,'FontName','Times','FontSize',16);
%     legend(['GS1: $ \ \mu_1 =' num2str(round(ad_1,1)) '\mathrm{\mu m}, \ \sigma_1 =' num2str(round(std_1,1)) '\mathrm{\mu m}$' ],...
%         ['GS2: $ \ \mu_2 =' num2str(round(ad_2,1)) '\mathrm{\mu m}, \ \sigma_2 =' num2str(round(std_2,1)) '\mathrm{\mu m}$' ],'location','northeast','Interpreter','latex','FontSize', 12);
%     % title(['$g1-' num2str(g1*1E6) '\ g2-' num2str(g2*1E6) '\ ' num2str(round(sample_width*1e3))...
%     %     '-' num2str(round(sample_depth_1*1e3)) '-' num2str(round(sample_depth_2*1e3)) '\ ' num2str(run_num) '$'],'Interpreter','latex');
% 
%     saveas(gcf, 'PDF_S1.png', 'png');
% 
%     % fwave = figure('Visible','on');
%     % histogram(grain_size_1*1e6,'binwidth',1,'Normalization', 'pdf'); hold on; %'FaceColor', [0, 0.2, 0.5]
%     % histogram(grain_size_2*1e6,'binwidth',1,'Normalization', 'pdf'); %'FaceColor',[0.5, 0.5, 0.5]
%     % xlabel('Grain Size ($\mu m$)','Interpreter','latex');
%     % ylabel('PDF','Interpreter','latex');
%     % box on;
%     % set(gca,'FontName','Times','FontSize',16);
%     % legend(['GS1: $ \ \mu =' num2str(round(ad_1*1E6,1)) '\mu m, \ \sigma =' num2str(round(std_1*1E6,1)) '\mu m$' ],...
%     %     ['GS2: $ \ \mu =' num2str(round(ad_2*1E6,1)) '\mu m, \ \sigma =' num2str(round(std_2*1E6,1)) '\mu m$' ],'location','northeast','Interpreter','latex');
%     % title(['$g1-' num2str(g1*1E6) '\ g2-' num2str(g2*1E6) '\ ' num2str(round(sample_width*1e3))...
%     %     '-' num2str(round(sample_depth_1*1e3)) '-' num2str(round(sample_depth_2*1e3)) '\ ' num2str(run_num) '$'],'Interpreter','latex');
%     %% Plot Image
% %     I2 = I1;
% % for i=1:max(max(I1))
% %     iI2 = randi([1,max(max(I1))*10]);
% %     I2(I2==i)=iI2;
% % end
% % figure;imagesc(I2);
%     %% Save data
%     save([sf_dr '\TPC' num2str(run_num) '.mat'], 'r_1', 'wr_1', 'r_2', 'wr_2','pixels_size');
%     save([sf_dr '\GRAINSIZE' num2str(run_num) '.mat'], 'grain_size_1', 'grain_size_2');
%     save([sf_dr '\MEANGRAINPROPERITY' num2str(run_num) '.mat'], 'aL_1', 'ad_1', 'aL_2', 'ad_2');
end
end
end