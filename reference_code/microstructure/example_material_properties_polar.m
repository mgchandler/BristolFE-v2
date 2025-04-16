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

%%%% euler methods, orientation angles have an uniform distribution
%rand('seed',sum(230*clock));
rand('seed',1e6);
euler_angle(:,1)=2*pi*rand(length(GRAIN_label),1);
euler_angle(:,2)=acos(2*rand(length(GRAIN_label),1)-1);
euler_angle(:,3)=2*pi*rand(length(GRAIN_label),1);

CA=fn_euler_c(C0,euler_angle(:,1),euler_angle(:,2),euler_angle(:,3));

for i=1:length(GRAIN_label)
%for i = 1:5
    % temp_material.material{i}.paramsType=anisotropic;
    temp_material.material{i}.paramValues=[CA{i}([1,7,8,13,14,15,19,20,21,22,25,26,27,28,29,31,32,33,34,35,36]) density];
    temp_material.matTypeRefs(GRAIN_label{i})=i;
end

% Convert Euler angles to crystallographic axes
% Assuming [phi1, Phi, phi2] convention for Euler angles
% Convert to [theta, phi] for stereographic projection
theta_h110 = atan2(sin(euler_angle(:,1)).*sin(euler_angle(:,2)), cos(euler_angle(:,1)));
phi_h110 = acos(cos(euler_angle(:,1)).*sin(euler_angle(:,2)));

theta_h111 = atan2(sin(euler_angle(:,3)).*cos(euler_angle(:,2)), sin(euler_angle(:,1)).*sin(euler_angle(:,2)) + cos(euler_angle(:,1)).*cos(euler_angle(:,2)).*cos(euler_angle(:,3)));
phi_h111 = acos(-sin(euler_angle(:,1)).*sin(euler_angle(:,2)).*cos(euler_angle(:,3)) + cos(euler_angle(:,1)).*cos(euler_angle(:,2)));

% Plot pole figure for h110 axis
figure;
subplot(1,2,1);
plot_pole_figure(theta_h110, phi_h110, 'h110 Pole Figure');

% Plot inverse pole figure for h110 axis
subplot(1,2,2);
plot_pole_figure(-theta_h110, phi_h110, 'Inverse h110 Pole Figure');

% Plot pole figure for h111 axis
figure;
subplot(1,2,1);
plot_pole_figure(theta_h111, phi_h111, 'h111 Pole Figure');

% Plot inverse pole figure for h111 axis
subplot(1,2,2);
plot_pole_figure(-theta_h111, phi_h111, 'Inverse h111 Pole Figure');

function plot_pole_figure(theta, phi, title_str)
    % Convert to Cartesian coordinates
    x = sin(phi) .* cos(theta);
    y = sin(phi) .* sin(theta);
    z = cos(phi);

    % Plot sphere
    [x_sphere, y_sphere, z_sphere] = sphere(100);
    surf(x_sphere, y_sphere, z_sphere, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    hold on;

    % Plot points
    plot3(x, y, z, 'b.', 'MarkerSize', 10);
    axis equal;
    view(3);
    title(title_str);
    xlabel('x');
    ylabel('y');
    zlabel('z');
end

