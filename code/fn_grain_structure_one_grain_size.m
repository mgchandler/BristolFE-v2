function [VX_1,VY_1,V_1,C_1]=fn_grain_structure_one_grain_size(mean_grain_size_1,...
    grain_size_spread_1,sample_width,sample_depth)


% 1
x_size=sample_width;
z_size_1=sample_depth;
x_start_1=-mean_grain_size_1*5:mean_grain_size_1:x_size+mean_grain_size_1*5; % creates 
z_start_1=-mean_grain_size_1*5:mean_grain_size_1:sample_depth+mean_grain_size_1*5;

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

% TEST
% figure;
% plot(grain_x_1,grain_z_1, 'o');  % Plot input points
% hold on;
% plot(VX_1,VY_1, '-'); 
% title('Voronoi Diagram');
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