function [matls] = fn_gen_grain_matls(elastic_matrix, rho, col, name, el_type,GRAIN_label)

C0=zeros(6,6);
    C0([1,7,8,13,14,15,19,20,21,22,25,26,27,28,29,31,32,33,34,35,36])=elastic_matrix;
    C0=C0+tril(C0',-1);
    
    %%%% euler methods, orientation angles have an uniform distribution
    % rand('seed',sum(230*clock));
    rand('seed',1e6);
    euler_angle(:,1)=2*pi*rand(length(GRAIN_label),1);
    euler_angle(:,2)=acos(2*rand(length(GRAIN_label),1)-1);
    euler_angle(:,3)=2*pi*rand(length(GRAIN_label),1);
    
    CA=fn_euler_c(C0,euler_angle(:,1),euler_angle(:,2),euler_angle(:,3));
    
    for i=1:length(GRAIN_label)
         %CPE3 must be the element type for a solid
        temp_material.material{i}.paramValues = CA{i};
        temp_material.matTypeRefs(GRAIN_label{i})=i;
    end

    for i = 1:length(temp_material.matTypeRefs)
        material_index = temp_material.matTypeRefs(i);
        temp_material.vals{i}.temp_D = temp_material.material{material_index};
        temp_material.vals{i}.rho = [rho];
        temp_material.vals{i}.col = [col];
        temp_material.vals{i}.name = [name];
        temp_material.vals{i}.el_type = [el_type];
        temp_material.vals{i}.D = temp_material.vals{i}.temp_D.paramValues;
    end

%--------------------------------------------------------------------------
%DEFINE THE PROBLEM

matls = temp_material.vals;

end