function [grains]=fn_label_grains_two_grain_size_using_VC(ele_cen,sample_depth_2,V_1,C_1,V_2,C_2)

tic

%
ele_cen_x=ele_cen(:,1);
ele_cen_z=ele_cen(:,2);
ele_ind_1=find(ele_cen_z>=sample_depth_2);
ele_ind_2=find(ele_cen_z<sample_depth_2);

ele_cen_x_1=ele_cen_x(ele_ind_1);
ele_cen_z_1=ele_cen_z(ele_ind_1);
ele_cen_x_2=ele_cen_x(ele_ind_2);
ele_cen_z_2=ele_cen_z(ele_ind_2);


in_ind_1=ele_ind_1;
V_x_1=V_1(:,1);
V_z_1=V_1(:,2);
grains_1=cell(length(C_1),1);
for i=1:length(C_1)
%     fn_show_progress(i,length(C_1),'Assigning grains',1,1);
    if sum(V_1(C_1{i},:)==inf)>0
    elseif sum(V_1(C_1{i},:)==-inf)>0
    else
        in=inpolygon(ele_cen_x_1,ele_cen_z_1,V_x_1(C_1{i}),V_z_1(C_1{i}));% in=1 for elementcentre in poly made up with V(C{i})
        grains_1{i}=in_ind_1(in);
    end
end

in_ind_2=ele_ind_2;
V_x_2=V_2(:,1);
V_z_2=V_2(:,2);
grains_2=cell(length(C_2),1);
for i=1:length(C_2)
%     fn_show_progress(i,length(C_2),'Assigning grains',1,1);
    if sum(V_2(C_2{i},:)==inf)>0
    elseif sum(V_2(C_2{i},:)==-inf)>0
    else
        in=inpolygon(ele_cen_x_2,ele_cen_z_2,V_x_2(C_2{i}),V_z_2(C_2{i}));% in=1 for elementcentre in poly made up with V(C{i})
        grains_2{i}=in_ind_2(in);
    end
end

grains=[grains_1;grains_2];

% Assuming your cell array is called 'yourCellArray'
notEmptyCells1 = ~cellfun(@isempty, grains_1); % Logical array indicating non-empty cells
numNotEmptyCells1 = sum(notEmptyCells1); % Count of non-empty cells

% Assuming your cell array is called 'yourCellArray'
notEmptyCells2 = ~cellfun(@isempty, grains_2); % Logical array indicating non-empty cells
numNotEmptyCells2 = sum(notEmptyCells2); % Count of non-empty cells
% mat2poly
% write_to_file

toc

end