close all
clear all

% crosstalk window 
xtalk_start=-3e-5;
xtalk_end=3e-5;
% setting up a window of time within which to look for cross-talk


load('Example_pristine_and_damage_ANTH_MLE.mat') % damanged and undamanged signal data
load('time.mat') %time values corresponding to above
max_t=0.5e-3; % max time value to 0.0005 secs

in_t=find(time<=max_t); 
time=time(in_t);% used to filter out time beyond a certain threshold (max time)


not_xtalk=find(time<xtalk_start | time>xtalk_end); %indices outside of crosstalk window

t_off=length(time)-length(not_xtalk)+1;

shift = -time(1);
time = time + shift;

undam=BAE_12_pristine(in_t);
dam=BAE_12_damage_pos1(in_t);

sub_sig=dam-undam; %difference between damaged and undamaged
processed=fn_rlt1(abs(sub_sig(not_xtalk)));

analytic_signal = hilbert(sub_sig);
hilbert_env = abs(analytic_signal);
% figure
% plot(time,real([dam undam]))

figure;
%plot(time,real(sub_sig)./max(real(sub_sig)))
hold on
plot(time(t_off:end)*10e2,real(sub_sig(not_xtalk))./max(real(sub_sig(not_xtalk))), 'Color', [0, 0.25, 0.5]);
plot(time(t_off:end)*10e2, real(hilbert_env(not_xtalk))./max(real(hilbert_env(not_xtalk))), 'Color', [0.4, 0.4, 0.4],'Linewidth',1.5 );
%plot(time(t_off:end),processed./max(processed));
xlabel('time (ms)', 'Interpreter', 'latex');
ylabel('Normalised Amplitude (-)', 'Interpreter', 'latex');
    
set(gca,'FontName','Times','FontSize',14);

% % Set the figure size such that the x-axis is twice as long as the y-axis
fig = gcf; % Get current figure handle
% current_position = fig.Position; % Get current figure position
% fig.Position = [current_position(1:2), current_position(3)*1.7, current_position(4)]; % Adjust figure position

% Save the figure
saveas(fig, 'signal.png'); % You can change the file format if needed (e.g., 'sample_figure.pdf', 'sample_figure.jpg')


figure;
plot(time(t_off:end)*10e2,processed./max(processed),  'Color', [0, 0.25, 0.5],'Linewidth',1.5);
xlabel('time (ms)', 'Interpreter', 'latex');
ylabel('Normalised filter result (-)', 'Interpreter', 'latex');
set(gca,'FontName','Times','FontSize',14);

fig2 = gcf; % Get current figure handle
% current_position = fig2.Position; % Get current figure position
% fig2.Position = [current_position(1:2), current_position(3)*1.7, current_position(4)]; % Adjust figure position

% Save the figure
saveas(fig2, 'rmle.png'); % You can change the file format if needed (e.g., 'sample_figure.pdf', 'sample_figure.jpg')


figure;
plot(time(t_off:end)*10e2,real(sub_sig(not_xtalk))./max(real(sub_sig(not_xtalk))));
hold on
plot(time(t_off:end)*10e2,processed./max(processed),'Linewidth',1.5);
xlabel('Time (ms)', 'Interpreter', 'latex');
ylabel('Normalised Amplitude (-)', 'Interpreter', 'latex');
set(gca,'FontName','Times','FontSize',14);
legend('Standardised Waveform', 'RMLE filter result', 'Interpreter', 'latex',  'FontSize', 14)
ylim([-1,1])

fig3 = gcf; % Get current figure handle
% current_position = fig3.Position; % Get current figure position
% fig3.Position = [current_position(1:2), current_position(3)*2, current_position(4)]; % Adjust figure position

% Save the figure
saveas(fig3, 'rmle_small.png'); % You can change the file format if needed (e.g., 'sample_figure.pdf', 'sample_figure.jpg')



%plot(time(t_off:end),processed);