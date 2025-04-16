% Sample data (replace with your actual data)
els_per_grain = [1.6666666, 2.727272, 3.3333333, 4.16666666, 5.1724, 6.5217, 7.5];  
v_error = [0.349, 0.281, 0.15, 0.131, 0.138, 0.105, 0];  
% a_error = [-11.672, -9.159, -8.261, -3.455, -2.121, -0.211, 0];
a_error = [-1.1672, -0.9159, -0.8261, -0.3455, -0.2121, -0.121, 0];




% Plotting
figure;
plot(els_per_grain, v_error, '.-','Color', [0, 0.25, 0.5], 'LineWidth', 2, 'MarkerSize', 25);
ylabel('$(V-V_{c})$ $V_{c}^{-1}$ (\%)', 'Interpreter', 'latex');
xlabel('$EG_{2}^{-1}$ (elements per grain)', 'Interpreter', 'latex');
set(gca,'FontName','Times','FontSize',12);

figure;
plot(els_per_grain, a_error, '.-', 'Color', [0, 0.25, 0.5], 'LineWidth', 2, 'MarkerSize', 25);
ylabel('$(A - A_{c}) A_{c}^{-1}$ (\%)', 'Interpreter', 'latex');
xlabel('$EG_{2}^{-1}$ (elements per grain)', 'Interpreter', 'latex');
set(gca,'FontName','Times','FontSize',12);

figure;
plot(els_per_grain, a_error, '.-', 'LineWidth', 2, 'MarkerSize', 25);
hold on
plot(els_per_grain, v_error, '.-', 'LineWidth', 2, 'MarkerSize', 25);

 legend(['$(A - A_{c}) A_{c}^{-1}$'],...
        ['$(t-t_{c})$ $t_{c}^{-1}$'],'location','southeast', 'Interpreter', 'latex',  'FontSize', 14);

ylabel('Error (\%)', 'Interpreter', 'latex');
xlabel('$EG_{1}^{-1}$ (elements per grain dimension)', 'Interpreter', 'latex');
set(gca,'FontName','Times','FontSize',14);

filename = 'Mesh_Convergence_3.png';
print(filename, '-dpng', '-r300');