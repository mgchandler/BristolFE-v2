function [matls, dis] = fn_gen_grain_matls_2(elastic_matrix, rho, col, name, el_type,GRAIN_label, set_col, el_size)

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

    % Define the colormap
    % Determine the number of unique values in temp_material.matTypeRefs
    unique_values = unique(temp_material.matTypeRefs);
    num_unique_values = numel(unique_values);

    unique_values = unique(temp_material.matTypeRefs);
    occurrences = histc(temp_material.matTypeRefs, unique_values);
    
    % Create a structure to store unique values and their counts
    unique_counts = struct();
    
    el_area = (el_size)^2 * sqrt(3) / 4; %area of one triangle element
    % Store unique values and their counts in the structure
    for i = 1:numel(unique_values)
        unique_counts(i).value = unique_values(i);
        unique_counts(i).count = occurrences(i);
        %unique_counts(i).grain_area = unique_counts(i).count * el_area;
    end

    % Calculate unique counts of occurrences
    unique_counts_of_occurrences = unique(occurrences);
    num_unique_counts = numel(unique_counts_of_occurrences);
    
    % Initialize structure to store unique counts and their corresponding number of instances
    unique_counts_data = struct('unique_count', {}, 'num_instances', {});
    
    % Loop through unique counts and calculate the number of instances
    for i = 1:num_unique_counts
        unique_count_value = unique_counts_of_occurrences(i);
        num_instances = sum(occurrences == unique_count_value);
        
        % Store data in structure
        unique_counts_data(i).unique_count = unique_count_value;
        unique_counts_data(i).num_instances = num_instances;
        unique_counts_data(i).grain_area = unique_counts_data(i).unique_count * el_area;
    end
    

    % Define the desired hue range (e.g., from 0.2 to 0.6)
    % hue_min = 0.1; %change when wanted
    % hue_max = 0.9; % change when wanted
    hue_min = 0.075; % Change when wanted
    hue_max = 0.925;
    
    % Calculate the number of hues needed within the specified range
    %num_hues = ceil((hue_max - hue_min) * num_unique_values);
    
    num_hues = num_unique_values;
    % Generate random hues within the specified range
    random_hues = linspace(hue_min, hue_max, num_hues);
    
    % Set constant values for saturation and value
    saturation = 0.75;
    value = 0.85;

    % saturation = 1;
    % value = 0.9;
    
    % Create the colormap using the generated hues, constant saturation, and value
    cmap = hsv2rgb([random_hues', saturation * ones(num_hues, 1), value * ones(num_hues, 1)]);
    
    % Randomly shuffle the rows of cmap
    cmap = cmap(randperm(size(cmap, 1)), :);
    % unique_values = unique(temp_material.matTypeRefs);
    % num_unique_values = numel(unique_values);
    % 
    % % Create a colormap with a size that matches the number of unique values
    % cmap = hsv(num_unique_values);
    

    for i = 1:length(temp_material.matTypeRefs)
        material_index = temp_material.matTypeRefs(i);
        idx = find(unique_values == material_index);
        %temp_material.vals{i}.col = cmap(idx, :);
        temp_material.vals{i}.MatRef = material_index;
        temp_material.vals{i}.temp_D = temp_material.material{material_index};
        temp_material.vals{i}.rho = [rho];
        %temp_material.vals{i}.col = [col];
        temp_material.vals{i}.name = [name];
        temp_material.vals{i}.el_type = [el_type];
        temp_material.vals{i}.D = temp_material.vals{i}.temp_D.paramValues;
        if exist('set_col', 'var') && set_col == 0
            temp_material.vals{i}.col = [col];
        else
            temp_material.vals{i}.col = cmap(idx, :);
        end

     

    end


%--------------------------------------------------------------------------
%DEFINE THE PROBLEM

dis = unique_counts_data;
matls = temp_material.vals;

end