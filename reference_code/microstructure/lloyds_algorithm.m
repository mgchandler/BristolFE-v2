% Generate random points
% Define the domain
x_min = 0;  % Minimum x value
x_max = 10; % Maximum x value
y_min = 0;  % Minimum y value
y_max = 10; % Maximum y value

% Generate random points within the specified domain
num_points = 1000; % Number of random points
x = (x_max - x_min) * rand(num_points, 1) + x_min;
y = (y_max - y_min) * rand(num_points, 1) + y_min;

% Compute Voronoi diagram before Lloyd's algorithm
[vx_before, vy_before] = voronoi(x, y);
cell_areas_before = compute_cell_areas(x, y, vx_before, vy_before);
std_dev_area_before = std(cell_areas_before);
cv_area_before = std_dev_area_before / mean(cell_areas_before);
regular_polygon_areas_before = (3 * sqrt(3) / 2) * (sqrt(cell_areas_before));
regular_polygon_ratio_before = cell_areas_before ./ regular_polygon_areas_before;

% Number of iterations for Lloyd's algorithm
num_iterations = 1000;
% Perform Lloyd's algorithm
for iter = 1:num_iterations
    % Compute Voronoi diagram
    [vx, vy] = voronoi(x, y);
    
    % Compute centroids of Voronoi cells
    centroids_x = zeros(numel(x), 1);
    centroids_y = zeros(numel(y), 1);
    for i = 1:numel(x)
        % Find indices of vertices for i-th cell
        cell_vertices = vy(:, i) ~= 0;
        
        % Compute centroid of the i-th cell
        centroids_x(i) = mean(vx(cell_vertices, i));
        centroids_y(i) = mean(vy(cell_vertices, i));
    end
    
    % Update input points to centroids
    x = centroids_x;
    y = centroids_y;
end


% Compute Voronoi diagram after Lloyd's algorithm
[vx_after, vy_after] = voronoi(x, y);
cell_areas_after = compute_cell_areas(x, y, vx_after, vy_after);
std_dev_area_after = std(cell_areas_after);
cv_area_after = std_dev_area_after / mean(cell_areas_after);
regular_polygon_areas_after = (3 * sqrt(3) / 2) * (sqrt(cell_areas_after));
regular_polygon_ratio_after = cell_areas_after ./ regular_polygon_areas_after;

% Display changes in uniformity metrics
disp('Changes in Uniformity Metrics after Lloyd''s Algorithm:');
disp(['Standard Deviation of Cell Areas: ', num2str(std_dev_area_before), ' -> ', num2str(std_dev_area_after)]);
disp(['Coefficient of Variation of Cell Areas: ', num2str(cv_area_before), ' -> ', num2str(cv_area_after)]);
disp(['Average Regular Polygon Ratio: ', num2str(mean(regular_polygon_ratio_before)), ' -> ', num2str(mean(regular_polygon_ratio_after))]);

% Visualize Voronoi diagram after Lloyd's algorithm
figure;
plot(x, y, 'ro');
hold on;
for i = 1:numel(x)
    % Find indices of vertices for i-th cell
    cell_vertices = vy_after(:, i) ~= 0;
    
    % Plot vertices of i-th cell
    plot(vx_after(cell_vertices, i), vy_after(cell_vertices, i), 'b-');
end
title('Voronoi Diagram after Lloyd''s Algorithm');
xlabel('X');
ylabel('Y');
axis equal;


% Function to compute cell areas
function cell_areas = compute_cell_areas(x, y, vx, vy)
    num_points = numel(x);
    cell_areas = zeros(num_points, 1);
    for i = 1:num_points
        % Find the indices of the vertices for the i-th cell
        cell_vertices = vy(:, i) ~= 0;
        
        % Compute the area of the i-th cell
        cell_areas(i) = polyarea(vx(cell_vertices, i), vy(cell_vertices, i));
    end
end


