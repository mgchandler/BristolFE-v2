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