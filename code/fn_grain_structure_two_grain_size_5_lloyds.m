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

% % Initial Voronoi diagram
[VX_1_before, VY_1_before] = voronoi(grain_x_1, grain_z_1);

figure;
plot(grain_x_1,grain_z_1, 'b.');  % Plot input points
hold on;
plot(VX_1_before,VY_1_before, 'k-'); 
title('Voronoi Before');
axis equal;
ylim([1.1e-3, 1.3e-3]);
xlim([1e-4, 2.8e-4]);
% hold off;

% Number of iterations for Lloyd's algorithm
num_iterations = 8;

% Perform Lloyd's algorithm
for iter = 1:num_iterations
    % Compute Voronoi diagram based on current generating points
    [VX_1, VY_1] = voronoi(grain_x_1, grain_z_1);
    
    % Compute centroids of Voronoi cells
    centroids_x = zeros(length(grain_x_1), 1);
    centroids_z = zeros(length(grain_z_1), 1);
    for i = 1:length(grain_x_1)
        % Find the vertices of the current Voronoi cell
        cell_vertices_idx = C_1{i};
        cell_vertices_x = V_1(cell_vertices_idx, 1);
        cell_vertices_z = V_1(cell_vertices_idx, 2);
        
         % Check for NaN or infinite values in cell vertices
        if any(isnan(cell_vertices_x)) || any(isnan(cell_vertices_z)) || ...
           any(isinf(cell_vertices_x)) || any(isinf(cell_vertices_z))
            % Replace infinite centroid values with generating points
            warning('Infinite values found in Voronoi cell vertices. Using generating points instead.');
            centroids_x(i) = grain_x_1(i);
            centroids_z(i) = grain_z_1(i);
        else
            % Compute centroid of the current Voronoi cell
            centroids_x(i) = mean(cell_vertices_x);
            centroids_z(i) = mean(cell_vertices_z);
       
        end
    end
end
% Update Voronoi diagram
[VX_1_after, VY_1_after] = voronoi(centroids_x, centroids_z);

figure;
plot(centroids_x,centroids_z, 'b.');  % Plot input points
hold on;
plot(VX_1_after,VY_1_after, 'k-'); 
title('Voronoi After');
axis equal;
ylim([1.15e-3, 1.3e-3]);
xlim([-0.5e-4, 0.8e-4]);



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
% figure;
% plot(grain_x_2,grain_z_2, 'r.');  % Plot input points
% hold on;
% plot(VX_2,VY_2, 'k-'); 
% %title('Voronoi Diagram');
% axis equal;
% 
% xlim([-6e-5, 6e-5]);
% ylim([1.2e-4, 3.2e-4]);
% hold off;
% 
% 
% % xticks([-6e-5, -3e-5, 0, 3e-5, 6e-5])
% % xticklabels({'0', '30', '60', '90', '120'})
% 
% xticks([-6e-5, -2e-5, 2e-5, 6e-5])
% xticklabels({'0', '40', '80', '120'})
% % Adjust y-axis tick labels
% yticks([1.2e-4, 1.6e-4, 2.0e-4, 2.4e-4, 2.8e-4, 3.2e-4]);
% %yticklabels({'0', '0.4e-4', '0.8e-4', '1.2e-4', '1.6e-4', '2.0e-4'});
% yticklabels({'0', '40', '80', '120', '160', '200'});
% set(gca,'FontName','Times','FontSize',16);
% xlabel('x-coordinate ($\mu m$)', 'Interpreter', 'latex');
% ylabel('y-coordinate ($\mu m$)', 'Interpreter', 'latex');
% 
% filename = 'G2_30_S5_V3.png';
% print(filename, '-dpng', '-r300'); % Specify the resolution (300 dpi in this case)
% 
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