function [mod, GRAIN_label] = fn_gen_grain_structure_and_mesh(G1,G2,RUN_NUM,elementtype,el_size,sample_width,sample_depth,sample_depth_1,sample_depth_2,mat)

for ig2=1:length(G2)
for ig1=1:length(G1)
g2 = G2(ig2); 
g1 = G1(ig1); 
for run_num=1:RUN_NUM(ig1)
    % Save outputs into folder for investigation
    % close all; clear I0;
    % filename = ['g1_' num2str(round(g1*1e6)) '_g2_' num2str(round(g2*1e6)) ...
    %     '_' num2str(round(sample_width*1e3)) '-' num2str(round(sample_depth_1*1e3)) '-' num2str(round(sample_depth_2*1e3))];
    % sf_dr = [sf_dr0 '\OUTPUT-grainstructure\' filename '\'];
    % mkdir(sf_dr);

    % correlated parameters: defines the anistropy of the grain spread. 
    grain_size_spread_1 = g1;
    grain_size_spread_2 = g2; 

    % Function to build voronoi structure
    [VX_1,VY_1,V_1,C_1,VX_2,VY_2,V_2,C_2]=fn_grain_structure_one_grain_size...
        (g1,grain_size_spread_1,g2,grain_size_spread_2,sample_width,sample_depth,sample_depth_1,sample_depth_2);
    % save([sf_dr  'grain_structure' num2str(run_num) '.mat'],...
    %         'VX_1','VY_1','V_1','C_1','VX_2','VY_2','V_2','C_2'); 

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
     [GRAIN_label] = fn_label_grains_two_grain_size_using_VC(ele_cen,sample_depth_2,V_1,C_1,V_2,C_2);
    
    % save([sf_dr 'grain_label' num2str(run_num) '_' elementtype '_size' num2str(round(el_size*1e9)) '.mat'],...
    %     'GRAIN_label');  
end
end
end