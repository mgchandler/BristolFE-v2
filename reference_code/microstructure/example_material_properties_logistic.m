clear all;
close all;
clc;

elastic_matrix = [2.0460 1.3770 2.0460 1.3770 1.3770 2.0460 0  0  0 1.2620  0       ...
    0 0  0  1.2620    0    0    0    0   0    1.2620]*1.0e+11; % stainless_steel
density = 7870; % information from edupack

number_of_grains = 1000;
for ii = 1 : number_of_grains
    GRAIN_label{ii} = round((rand(5,1)+1)*100);  %% need to be replaced by properly assigning finite elements to each grains.
end

C0=zeros(6,6);
C0([1,7,8,13,14,15,19,20,21,22,25,26,27,28,29,31,32,33,34,35,36])=elastic_matrix;
C0=C0+tril(C0',-1); % allocates values from elastic_matrix into C0 constructing the 6x6 matrix

% Set parameters for the logistic distribution of rotation angles
mu_rotation = 180; % Mean (midpoint of the range)
beta_rotation = 30; % Scale parameter (controls spread)

% Generate random values following a logistic distribution for rotation angles
logistic_values_rotation = mu_rotation + beta_rotation * log(rand(length(GRAIN_label), 1));

% Generate random Euler angles
rand('seed',1e6);
euler_angle(:,1) = 2*pi*rand(length(GRAIN_label),1);
euler_angle(:,2) = acos(2*rand(length(GRAIN_label),1)-1);
euler_angle(:,3) = 2*pi*rand(length(GRAIN_label),1);

% Calculate rotation angles from Euler angles
rotation_angles = sqrt(euler_angle(:,1).^2 + euler_angle(:,2).^2 + euler_angle(:,3).^2);

% Calculate the scaling factor for Euler angles to match the desired logistic distribution of rotation angles
scale_factor = logistic_values_rotation ./ rotation_angles;

% Scale the Euler angles
euler_angle_scaled = euler_angle .* scale_factor;

CA=fn_euler_c(C0,euler_angle(:,1),euler_angle(:,2),euler_angle(:,3));

for i=1:length(GRAIN_label)
%for i = 1:5
    % temp_material.material{i}.paramsType=anisotropic;
    temp_material.material{i}.paramValues=[CA{i}([1,7,8,13,14,15,19,20,21,22,25,26,27,28,29,31,32,33,34,35,36]) density];
    temp_material.matTypeRefs(GRAIN_label{i})=i;
end

euler_angles_deg = rad2deg(euler_angle_scaled);
rotation_angles_2 = sqrt(euler_angle_scaled(:,1).^2 + euler_angle_scaled(:,2).^2 + euler_angle_scaled(:,3).^2);
%rotation_angles_deg = rad2deg(rotation_angles_2);

% Constrain rotation angles within the range of 0 to 180 degrees
rotation_angles_deg_constrained = mod(rotation_angles_2, 180);

% Plot the normalized probability distribution of the rotation angles
figure;
histogram(rotation_angles_2, 'Normalization', 'probability');
xlabel('Rotation Angle (degrees)');