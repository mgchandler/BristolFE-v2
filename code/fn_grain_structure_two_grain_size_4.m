function [VX_1,VY_1,V_1,C_1,VX_2,VY_2,V_2,C_2]=fn_grain_structure_two_grain_size(mean_grain_size_1,...
    grain_size_spread_1,mean_grain_size_2,grain_size_spread_2,sample_width,sample_depth,sample_depth_1,sample_depth_2)


% 1
x_size=sample_width;
z_size_1=sample_depth_1;
x_start_1=-mean_grain_size_1*5:mean_grain_size_1:x_size+mean_grain_size_1*5; 
z_start_1=sample_depth_2-mean_grain_size_1*5:mean_grain_size_1:sample_depth+mean_grain_size_1*5;

grain_x_1=zeros(length(x_start_1)*length(z_start_1),1);
grain_z_1=grain_x_1;

ii=1;
for i_x = x_start_1
    for i_z = z_start_1
        grain_x_1(ii)= i_x + grain_size_spread_1*randn(1);
        grain_z_1(ii)= i_z + grain_size_spread_1*randn(1);
        ii=ii+1;
    end
end
[VX_1,VY_1]=voronoi(grain_x_1,grain_z_1);  % the finite vertices of the Voronoi edges (a start and end points in vertices) in vx and vy, Delaunay triangulation,  Voronoi tesselation
[V_1,C_1]=voronoin([grain_x_1,grain_z_1]); % Voronoi vertices V and the Voronoi cells C of the Voronoi diagram of X.
% 
% Define the hue range and number of unique values
hue_min = 0.075; % Change when wanted
hue_max = 1; % Change when wanted
num_unique_values = length(VX_1);

% Calculate the number of hues needed within the specified range
num_hues = num_unique_values;

% Generate random hues within the specified range
random_hues = linspace(hue_min, hue_max, num_hues);

% Set constant values for saturation and value
saturation = 0.75;
value = 0.85;

% Create the colormap using the generated hues, constant saturation, and value
cmap = hsv2rgb([random_hues', saturation * ones(num_hues, 1), value * ones(num_hues, 1)]);

% Randomly shuffle the rows of cmap
cmap = cmap(randperm(size(cmap, 1)), :);


angles = rand(size(grain_x_1)) * 2 * pi;

% Define arrow length
arrow_length = 1e-5; % Adjust this value as needed



%TEST


% Generate random Euler angles for each cluster
num_clusters = length(C_1);
random_euler_angles = rand(num_clusters, 3) * 360; % Random Euler angles in degrees

figure;
plot(grain_x_1,grain_z_1, 'r.');  % Plot input points
hold on;
% for i = 1:length(C_1)
%     if all(C_1{i} ~= 1) && all(C_1{i} ~= Inf)
%         % Plot patch with color from colormap
%         patch(V_1(C_1{i}, 1), V_1(C_1{i}, 2), cmap(i,:), 'FaceAlpha', 0.9);
%     end
% end
% 
% % Ensure ticks are increasing and unique
% ticks = unique(linspace(0, 360, 5));
% ticks = [0, 45, 90, 135, 180, 225, 270, 315, 360];
% % Define the restricted colormap
% hue_min = 0.075;
% hue_max = 0.65;
% restricted_cmap = hsv(ceil(size(cmap, 1) * (hue_max - hue_min)));
% restricted_cmap = restricted_cmap(floor(size(cmap, 1) * hue_min):end, :);
% 
% % Create colorbar legend for Euler angles
% colormap(restricted_cmap); % Set the colormap to restricted_cmap
% c = colorbar; % Create the colorbar
% 
% % Set colorbar ticks and labels based on Euler angles
% c.Ticks = linspace(0, 1, length(ticks)); % Set colorbar ticks positions
% c.TickLabels = arrayfun(@(x) sprintf('%.0f°', x), ticks, 'UniformOutput', false); % Set colorbar tick labels
% c.Label.String = 'Euler Angles (degrees)'; % Set colorbar label

% for i = 1:length(C_1)
%     if all(C_1{i} ~= 1) && all(C_1{i} ~= Inf)
%         patch(V_1(C_1{i}, 1), V_1(C_1{i}, 2), rand(1,3), 'FaceAlpha', 0.3);
%         patch(V_1(C_1{i}, 1), V_1(C_1{i}, 2), cmap(i,:), 'FaceAlpha', 0.9); % Random color with transparency
%     end
% end

% % Plot arrows at each point
% for i = 1:numel(grain_x_1)
%     dx = arrow_length * cos(angles(i));
%     dz = arrow_length * sin(angles(i));
%     quiver(grain_x_1(i), grain_z_1(i), dx, dz, 0, 'color', 'w', 'MaxHeadSize', 1.5); % 'r' for red arrows, adjust color as needed
% end
% 
plot(VX_1,VY_1, 'k-'); 
%title('Voronoi Diagram');
axis equal;

%xlim([-6e-5, 6e-5]);
xlim([-6e-5, 6e-5]);
%ylim([2.2e-4, 4.2e-4]);
ylim([6.2e-4, 8.2e-4]);


hold off;

xticks([-6e-5, -2e-5, 2e-5, 6e-5])
xticklabels({'0', '0.04', '0.08', '0.12'})
% Adjust y-axis tick labels
%yticks([2.2e-4, 2.6e-4, 3.0e-4, 3.4e-4, 3.8e-4, 4.2e-4]);
yticks([6.2e-4, 6.6e-4, 7.0e-4, 7.4e-4, 7.8e-4, 8.2e-4]);
yticklabels({'0', '0.04', '0.08', '0.12', '0.16', '0.20'});
set(gca,'FontName','Times','FontSize',14);

xlabel('x-coordinate (mm)', 'Interpreter', 'latex');
ylabel('y-coordinate (mm)', 'Interpreter', 'latex');
% 
% % Save the figure as a PNG file with high resolution
filename = 'G1_15_S5_SeedVor_14.png';
print(filename, '-dpng', '-r300'); % Specify the resolution (300 dpi in this case)

% 2
x_size=sample_width;
z_size_2=sample_depth_2;
x_start_2=-mean_grain_size_2*5:mean_grain_size_2:x_size+mean_grain_size_2*5;
z_start_2=-mean_grain_size_2*5:mean_grain_size_2:sample_depth_2+mean_grain_size_2*5;

grain_x_2=zeros(length(x_start_2)*length(z_start_2),1);
grain_z_2=grain_x_2;

ii=1;
for i_x = x_start_2
    for i_z = z_start_2
        grain_x_2(ii)= i_x + grain_size_spread_2*randn(1);
        grain_z_2(ii)= i_z + grain_size_spread_2*randn(1);
        ii=ii+1;
    end
end
[VX_2,VY_2]=voronoi(grain_x_2,grain_z_2);  % the finite vertices of the Voronoi edges (a start and end points in vertices) in vx and vy, Delaunay triangulation,  Voronoi tesselation
[V_2,C_2]=voronoin([grain_x_2,grain_z_2]); % Voronoi vertices V and the Voronoi cells C of the Voronoi diagram of X.

%TEST
figure;
plot(grain_x_2,grain_z_2, 'r.');  % Plot input points
hold on;
plot(VX_2,VY_2, 'k-'); 
%title('Voronoi Diagram');
axis equal;

xlim([-6e-5, 6e-5]);
ylim([1.2e-4, 3.2e-4]);
hold off;


% xticks([-6e-5, -3e-5, 0, 3e-5, 6e-5])
% xticklabels({'0', '30', '60', '90', '120'})

xticks([-6e-5, -2e-5, 2e-5, 6e-5]);
xticklabels({'2.50', '2.54', '2.58', '2.62'});
% Adjust y-axis tick labels
yticks([1.2e-4, 1.6e-4, 2.0e-4, 2.4e-4, 2.8e-4, 3.2e-4]);
%yticklabels({'0', '0.4e-4', '0.8e-4', '1.2e-4', '1.6e-4', '2.0e-4'});
yticklabels({'0', '0.04', '0.08', '0.12', '0.16', '0.20'});
set(gca,'FontName','Times','FontSize',14);
xlabel('x-coordinate (mm)', 'Interpreter', 'latex');
ylabel('y-coordinate (mm)', 'Interpreter', 'latex');
filename = 'G2_30_S5_14.png';
print(filename, '-dpng', '-r300'); % Specify the resolution (300 dpi in this case)

% % 
% %TEST
% figure;
% plot(grain_x_1,grain_z_1, 'o');  % Plot input points
% hold on;
% plot(VX_1,VY_1, '-'); 
% hold on;
% plot(grain_x_2,grain_z_2, 'o');  % Plot input points
% hold on;
% plot(VX_2,VY_2, '-'); 
% title('Voronoi Diagram GS1');
% axis equal;
% hold off;


% grain_structure.VX_1 = VX_1;
% grain_structure.VY_1 = VY_1;
% grain_structure.V_1 = V_1;
% grain_structure.C_1 = C_1;
% grain_structure.VX_2 = VX_2;
% grain_structure.VY_2 = VY_2;
% grain_structure.V_2 = V_2;
% grain_structure.C_2 = C_2;


end