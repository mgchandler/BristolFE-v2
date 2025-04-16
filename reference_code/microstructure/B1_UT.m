%%% Generate random grain structures with two different grain sizes,
%%% and calculate the TPCF and grain size data.
%%% Then, plot TPC curve and grain size distribution.
%%%%%%%%%
% close all;
clear;
% G1=[ 28 26 24 22 20 15]*1e-6;
G1=[15]*1e-6;
G2=[30]*1e-6;
% RUN_NUM=[100,100,100,100,100,100];
RUN_NUM = [1];
sample_width = 3e-3;
sample_depth = 10e-3;
sample_depth_1 = 4e-3;
sample_depth_2 = sample_depth-sample_depth_1;
numrun = 1e6;
element_size = 5e-6;
TPC_r_step = 3e-6;
elementtype = 'CPE4';
sf_dr0 = [ pwd '\'];
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
    grain_size_spread_1 = g1/5;
    grain_size_spread_2 = g2/5;    
    %% create new set of random grains
    [VX_1,VY_1,V_1,C_1,VX_2,VY_2,V_2,C_2]=fn_grain_structure_two_grain_size...
        (g1,grain_size_spread_1,g2,grain_size_spread_2,sample_width,sample_depth,sample_depth_1,sample_depth_2);
    save([sf_dr  'grain_structure' num2str(run_num) '.mat'],...
            'VX_1','VY_1','V_1','C_1','VX_2','VY_2','V_2','C_2');
    %% generate mesh
    if isfile([sf_dr0 '\OUTPUT-grainstructure\MESH\' elementtype  '_size' num2str(round(element_size*1e9)) '_' num2str(sample_width*1e6) '-' num2str(sample_depth*1e6) '.mat']);
        load([sf_dr0 '\OUTPUT-grainstructure\MESH\' elementtype  '_size' num2str(round(element_size*1e9)) '_' num2str(sample_width*1e6) '-' num2str(sample_depth*1e6) '.mat']);
    else
        if ismember(elementtype,{'CPE4','CPS4','CPE4R','CPS4R'})
            [nodes,elements,columns,rows] = fn_make_quad_mesh(element_size,sample_width,sample_depth);
            ele_cen=[(nodes(elements(:,1),1)+nodes(elements(:,2),1)+nodes(elements(:,3),1)+nodes(elements(:,4),1))/4,...
                (nodes(elements(:,1),2)+nodes(elements(:,2),2)+nodes(elements(:,3),2)+nodes(elements(:,4),2))/4];
        end
        mkdir([sf_dr0 '\OUTPUT-grainstructure\MESH\']);
        save([sf_dr0 '\OUTPUT-grainstructure\MESH\' elementtype '_size' num2str(round(element_size*1e9)) '_' num2str(sample_width*1e6) '-' num2str(sample_depth*1e6) '.mat'],...
            'nodes','elements','ele_cen','columns','rows');
    end
    %% Assigns the mesh elements to each grain
    GRAIN_label=fn_label_grains_two_grain_size_using_VC(ele_cen,sample_depth_2,V_1,C_1,V_2,C_2);
    save([sf_dr 'grain_label' num2str(run_num) '_' elementtype '_size' num2str(round(element_size*1e9)) '.mat'],...
        'GRAIN_label');   
    %% convert image data 
    rows_1 = round(rows*sample_depth_1/sample_depth_2);
    for i=1:length(GRAIN_label)
        I0(GRAIN_label{i})=i;
    end
    I1=reshape(I0,columns,rows);
    pixels_size = sample_depth/rows;%% 1 pixels 
    %% Measure TPC points for GRAIN SIZE 1
    I2_1 = I1(:,rows_1+1:end);
    R_1 = [TPC_r_step:TPC_r_step:g1*2]; 
    r_1 = unique(round(R_1/pixels_size)); % unit: pixels   
    [r_1, wr_1] = fn_calcTPCImage(I2_1, r_1, numrun);
    r_1  = [0,r_1];
    wr_1 = [1,wr_1];
    %% Measure TPC points for GRAIN SIZE 2
    I2_2 = I1(:,1:rows_1);
    R_2 = [TPC_r_step:TPC_r_step:g2*2];
    r_2 = unique(round(R_2/pixels_size)); % unit: pixels
    [r_2, wr_2] = fn_calcTPCImage(I2_2, r_2, numrun);
    r_2  = [0,r_2];
    wr_2 = [1,wr_2];
    %% Calculate the correlation length from the TPC curve
    ind_1 = find(wr_1 >= 0.1);
    aL_1 = interp1(wr_1(ind_1), r_1(ind_1)*pixels_size, exp(-1));
    ind_2 = find(wr_2 >= 0.1);
    aL_2 = interp1(wr_2(ind_2), r_2(ind_2)*pixels_size, exp(-1));
    %% PlotTPC curve
    fwave = figure('Visible','on');
    plot(r_1*pixels_size*1e6, wr_1,'Linewidth',2); hold on;
    plot(r_2*pixels_size*1e6, wr_2,'Linewidth',2);
    xlabel('r (\mum)');
    ylabel('TPCF');
    box on;
    set(gca,'FontName','Times','FontSize',16);
    title(['g1-' num2str(g1*1E6) ' g2-' num2str(g2*1E6) ' ' num2str(round(sample_width*1e3))...
        '-' num2str(round(sample_depth_1*1e3)) '-' num2str(round(sample_depth_2*1e3)) ' ' num2str(run_num)]);
    legend(['g1-' num2str(g1*1E6) '  \muL1-' num2str(round(aL_1*1E6,1))],...
        ['g2-' num2str(g2*1E6) '  \muL2-' num2str(round(aL_2*1E6,1))],'location','best');
    print(fwave, [sf_dr '\TPC' num2str(run_num) '.emf'], '-dmeta');
    %% Measure grain size distribution for GRAIN SIZE 1
    [grain_size_1] = fn_calcGrainSizeImage(I2_1);
    grain_size_1 = grain_size_1*sqrt((sample_depth/rows)*(sample_width/columns));
    %% Measure grain size distribution for GRAIN SIZE 2
    [grain_size_2] = fn_calcGrainSizeImage(I2_2);
    grain_size_2 = grain_size_2*sqrt((sample_depth/rows)*(sample_width/columns));
    %% Calculate the mean grainsize 
    ad_1 = mean(grain_size_1);
    ad_2 = mean(grain_size_2);
    %% Plot grain size distribution
    fwave = figure('Visible','on');
    histogram(grain_size_1*1e6,'binwidth',1,'Normalization', 'pdf', 'FaceColor', [0.2, 0.6, 0.8]); hold on; 
    histogram(grain_size_2*1e6,'binwidth',1,'Normalization', 'pdf' ,'FaceColor',[0.8, 0.2, 0.2]);
    xlabel('Grain Size(\mum)');
    ylabel('PDF');
    box on;
    set(gca,'FontName','Times','FontSize',16);
    legend(['g1-' num2str(g1*1E6) '  \mud1-' num2str(round(ad_1*1E6,1))],...
        ['g2-' num2str(g2*1E6) '  \mud2-' num2str(round(ad_2*1E6,1))],'location','best');
    title(['g1-' num2str(g1*1E6) ' g2-' num2str(g2*1E6) ' ' num2str(round(sample_width*1e3))...
        '-' num2str(round(sample_depth_1*1e3)) '-' num2str(round(sample_depth_2*1e3)) ' ' num2str(run_num)]);
    print(fwave, [sf_dr '\PDF' num2str(run_num) '.emf'], '-dmeta');
    %% Save data
    save([sf_dr '\TPC' num2str(run_num) '.mat'], 'r_1', 'wr_1', 'r_2', 'wr_2','pixels_size');
    save([sf_dr '\GRAINSIZE' num2str(run_num) '.mat'], 'grain_size_1', 'grain_size_2');
    save([sf_dr '\MEANGRAINPROPERITY' num2str(run_num) '.mat'], 'aL_1', 'ad_1', 'aL_2', 'ad_2');
end
end
end