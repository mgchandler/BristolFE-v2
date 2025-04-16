% Load data from the text file
data = load('conv_2.txt');

% Extract columns
x = data(:, 1);  % Column 1 for x-axis
y1 = data(:, 2); % Column 2 for y-axis
y2 = data(:, 3); % Column 3 for y-axis

% Plot data points
scatter(x, y1, 'b', 'filled'); hold on; % Plot y1
scatter(x, y2, 'r', 'filled'); % Plot y2

% Add labels and legend
xlabel('Column 1');
ylabel('Values');


% Show plot
grid on; % Show grid

figure;
scatter(x, y1, 'filled'); hold on; % Plot y1
set(gca,'FontName','Times','FontSize',14);

% Add labels and legend
xlabel('Number of Realisations','Interpreter', 'latex');
ylabel('Boundary Location Estimate','Interpreter', 'latex');
ylim([0, 5]); % Set y-axis limits
yline(2.5, '--k')
xticks(0:4:20);
filename = ['I1_Comp.png'];
print(filename, '-dpng', '-r300');

figure;
scatter(x, y2, 'filled', 'MarkerFaceColor', [0.8500, 0.3250, 0.0980]); % Orange
set(gca,'FontName','Times','FontSize',14);

% Add labels and legend
xlabel('Number of Realisations','Interpreter', 'latex');
ylabel('Boundary Location Estimate','Interpreter', 'latex');
ylim([0, 5]); % Set y-axis limits
yline(2.5, '--k')
xticks(0:4:20);
filename = 'I5_Comp.png';
print(filename, '-dpng', '-r300');