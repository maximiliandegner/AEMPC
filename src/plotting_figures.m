%% figures_journal.m 
% This script plots the results for comparing AE-MPC to E-MPC (no
% adaptation), E-MPC (perfect parameter knowledge)
% The following figures will be created:
% - Figure(1): Closed-loop state-space evolution
% - Figure(13): Closed-loop evolution of state [x]_2 (--> \ell(x,u) = [0;-1;0]*x)
% - Figure(2): Comparison of performance of MPC schemes, *time-invariant* case
% - Figure(6): Comparison of performance of MPC schemes, *timevarying* case

lw = 1.5; fs = 15; 
lscon   = '--';         lstrue  = '-';
colcon  = 'k';          coltrue = [1, 0, 0, 0.6];
col1 = [0, 0, 1,0.5]; col4 = [0.85, 0, 0, 1]; col2 = [0.23,0.57,0.13];
grey_alpha = 0.5;
lw2 = 2; lw4 = 3;
ls2 = ':'; ls4 = '--';
TV = false;
theta_true = syst.theta_true;
theta_nom = syst.theta_nom;
run parameter_def.m

sw_time = [1500; 2500];      % Time of switching to theta_true time-varying
rampval = [0.00000025;-0.000004];
ramp_end = [rampval(1)*(sw_time(2)-sw_time(1)); rampval(2)*(sw_time(2)-sw_time(1))];

saving_fig = true;
include_title = false;
%%



%% figure(1):
% load TI parameter files, comparing AE- with two E-MPC schemes
% only plotting the state-behavior
%load("Results/055/LMI_designed/TV/AE_TV2--N_25-Tsim_1000-mu_150-h_0.025-LMS_1-TV_1.mat", "x_arr", "theta_arr", "Tsim")
load("Results/Simulations/AE_scaled--N_25-Tsim_2500-mu_15-h_0.025-LMS_1-TV_0.mat", "x_arr", "theta_arr", "Tsim", "beta")
%%%load("Results/Simulations/current_version/AE_scaled--N_25-Tsim_50-mu_15-h_0.025-LMS_1-TV_0.mat", "x_arr", "theta_arr", "Tsim")
x_arr_AE = x_arr;
theta_arr_AE = theta_arr;

% load("Results/055/LMI_designed/TV/E_wrongTV--N_25-Tsim_1000-mu_150-h_0.025-LMS_0-TV_1.mat", "x_arr", "theta_arr")
load("Results/Simulations/E_wrong--N_25-Tsim_2500-mu_15-h_0.025-LMS_0-TV_0.mat", "x_arr", "theta_arr")
%%%load("Results/Simulations/current_version/E_wrong--N_25-Tsim_50-mu_15-h_0.025-LMS_0-TV_0.mat", "x_arr", "theta_arr")
x_arr_Earts = x_arr;
theta_arr_Earts = theta_arr;

% load("Results/055/LMI_designed/TV/E_trueTV2--N_25-Tsim_4000-mu_150-h_0.025-LMS_0-TV_1.mat", "x_arr", "theta_arr")
load("Results/Simulations/E_true--N_25-Tsim_2500-mu_15-h_0.025-LMS_0-TV_0.mat", "x_arr", "theta_arr")
%%%load("Results/Simulations/current_version/E_true--N_25-Tsim_50-mu_15-h_0.025-LMS_0-TV_0.mat", "x_arr", "theta_arr")
x_arr_AEperf = x_arr;
theta_arr_AEperf = theta_arr;

% col1 = [0, 0, 1, 0.5]; col4 = [1, 0, 0, 0.5]; col2 = [0.5, 1, 0, 0.5];
lstrue = "-.";
% fs = fs+1;

percent = 0.087; % percentage of last and first time span relative to Tsim
% OLD coloring (brigther colors, all lines "-" style): 
% col1 = [0, 0, 1, 0.5]; col4 = [1, 0, 0, 0.5]; col2 = [0.5, 1, 0, 0.5];
% lw2 = lw; lw4 = lw; ls2 = "-"; ls4 = "-";
f = figure(1); t = tiledlayout(3,2, "TileSpacing","compact", "Padding","tight"); 
f.Position = [300 60 800 625];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1 = nexttile; %top left
plot(x_arr_AE(1,1:Tsim), 'Color', col1, "LineWidth", lw); hold on;
plot(x_arr_AEperf(1,1:end), 'Color', col2, "LineWidth", lw2, "LineStyle", ls2);
plot(x_arr_Earts(1,1:Tsim), 'Color', col4, "LineWidth", lw4, "LineStyle", ls4);
plot([200; 200], [0.03; 0.13], "LineWidth", lw, "LineStyle","--", 'Color',[0,0,0,grey_alpha]);
%xlabel("Time", "Interpreter", "latex", 'FontSize', fs);
ylabel("$$[x]_1$$", "Interpreter", "latex", 'FontSize', fs);
yticks([0.03, 0.05, 0.07, 0.09, 0.11, 0.13])
ylim([0.03, 0.11])
grid on;
set(gca,'fontsize',fs-2);
set(gca,'XTickLabel',[]);
set(gca, "TickLabelInterpreter", "latex");

ax2 = nexttile; % top right
plot(x_arr_AE(1,1:Tsim), 'Color', col1, "LineWidth", lw); hold on;
plot(x_arr_AEperf(1,1:end), 'Color', col2, "LineWidth", lw2, "LineStyle", ls2);
plot(x_arr_Earts(1,1:Tsim), 'Color', col4, "LineWidth", lw4, "LineStyle", ls4);
plot([200; 200], [0.03; 0.13], "LineWidth", lw, "LineStyle","--", 'Color',[0,0,0,grey_alpha]);
%xlabel("Time", "Interpreter", "latex", 'FontSize', fs);
%ylabel("$$[x]_2$$", "Interpreter", "latex", 'FontSize', fs);
yticks([0.03, 0.05, 0.07, 0.09, 0.11, 0.13])
ylim([0.03, 0.11])
grid on;
set(gca,'fontsize',fs-2);
set(gca,'XTickLabel',[]);
% set(gca,'YTickLabel',[]);
set(gca, 'YAxisLocation', 'right');
set(gca, "TickLabelInterpreter", "latex");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax3 = nexttile; %middle left
plot(x_arr_AE(2,1:Tsim), 'Color', col1, "LineWidth", lw); hold on;
plot(x_arr_AEperf(2,1:end), 'Color', col2, "LineWidth", lw2, "LineStyle", ls2);
plot(x_arr_Earts(2,1:Tsim), 'Color', col4, "LineWidth", lw4, "LineStyle", ls4);
plot([200; 200], [0.078; 0.092], "LineWidth", lw, "LineStyle","--", 'Color',[0,0,0,grey_alpha]);
%xlabel("Time", "Interpreter", "latex", 'FontSize', fs);
ylabel("$$[x]_3$$", "Interpreter", "latex", 'FontSize', fs);
% legend('AE-MPC closed-loop state evolution', 'E-MPC closed-loop state evolution', 'FontSize', 14, 'Orientation', 'horizontal', 'Interpreter', 'latex');
yticks([0.078, 0.082, 0.086, 0.090])
ylim([0.078, 0.092])
grid on;
set(gca,'fontsize',fs-2);
set(gca,'XTickLabel',[]);
set(gca, "TickLabelInterpreter", "latex");

ax4 = nexttile; %middle right
plot(x_arr_AE(2,1:Tsim), 'Color', col1, "LineWidth", lw); hold on;
plot(x_arr_AEperf(2,1:end), 'Color', col2, "LineWidth", lw2, "LineStyle", ls2);
plot(x_arr_Earts(2,1:Tsim), 'Color', col4, "LineWidth", lw4, "LineStyle", ls4);
plot([200; 200], [0.078; 0.092], "LineWidth", lw, "LineStyle","--", 'Color',[0,0,0,grey_alpha]);
%xlabel("Time", "Interpreter", "latex", 'FontSize', fs);
%ylabel("$$[x]_1$$", "Interpreter", "latex", 'FontSize', fs);
yticks([0.078, 0.082, 0.086, 0.090])
ylim([0.078, 0.092])
grid on;
set(gca,'fontsize',fs-2);
set(gca,'XTickLabel',[]);
% set(gca,'YTickLabel',[]);
set(gca, 'YAxisLocation', 'right');
set(gca, "TickLabelInterpreter", "latex");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax5 = nexttile; % bottom left
plot(x_arr_AE(3,1:Tsim), 'Color', col1, "LineWidth", lw); hold on;
plot(x_arr_AEperf(3,1:end), 'Color', col2, "LineWidth", lw2, "LineStyle", ls2);
plot(x_arr_Earts(3,1:Tsim), 'Color', col4, "LineWidth", lw4, "LineStyle", ls4);
plot([200; 200], [0.13; 0.2], "LineWidth", lw, "LineStyle","--", 'Color',[0,0,0,grey_alpha]);
xlabel("Time", "Interpreter", "latex", 'FontSize', fs);
ylabel("$$[x]_2$$", "Interpreter", "latex", 'FontSize', fs);
yticks(0.13:0.02:0.2)
ylim([0.13,0.2])
grid on;
set(gca,'fontsize',fs-2);
set(gca, "TickLabelInterpreter", "latex");

ax6 = nexttile; %bottom right
plot(x_arr_AE(3,1:Tsim), 'Color', col1, "LineWidth", lw); hold on;
plot(x_arr_AEperf(3,1:end), 'Color', col2, "LineWidth", lw2, "LineStyle", ls2);
plot(x_arr_Earts(3,1:Tsim), 'Color', col4, "LineWidth", lw4, "LineStyle", ls4);
plot([200; 200], [0.13; 0.2], "LineWidth", lw, "LineStyle","--", 'Color',[0,0,0,grey_alpha]);
xlabel("Time", "Interpreter", "latex", 'FontSize', fs);
%ylabel("$$[x]_3$$", "Interpreter", "latex", 'FontSize', fs);
leg = legend('AE-MPC', 'E-MPC with $$\theta^\ast$$', 'E-MPC with $$\hat{\theta}$$', 'FontSize', fs, 'Orientation', 'horizontal', 'Interpreter', 'latex');
leg.Layout.Tile = 'north';
yticks(0.13:0.02:0.2)
ylim([0.13,0.2])
grid on;
set(gca,'fontsize',fs-2);
% set(gca,'YTickLabel',[]);
set(gca, 'YAxisLocation', 'right');
set(gca, "TickLabelInterpreter", "latex");

linkaxes([ax1 ax3 ax5], 'x')
xlim(ax1, [0, percent*Tsim*1.05])
linkaxes([ax2 ax4 ax6], 'x')
xlim(ax2, [(1-percent)*Tsim*0.997, Tsim*1.003])

if include_title
    title(t, "Closed-Loop State Evolution", 'FontSize', fs+3, 'Interpreter', 'latex');
end

if saving_fig
    N = input('Horizon length:\n')
    beta = input('Value of beta:\n')
    LMS = 1;
    mu = input('Value of mu:\n')
    comment = input('Additional comments (may be empty):', 's')
    exportgraphics(gcf, ['Results/Figures/States_TI-param' num2str(N) '_' num2str(Tsim) '_beta_' num2str(beta) '_LMS_' num2str(LMS) '_mu_' num2str(mu) '_comment-' comment '.pdf'])
    exportgraphics(gcf, ['Results/Figures/States_TI-param' num2str(N) '_' num2str(Tsim) '_beta_' num2str(beta) '_LMS_' num2str(LMS) '_mu_' num2str(mu) '_comment-' comment '.eps'])
end


%% figure(13): State x2 at beginning and end of the simulation
f = figure(13); t = tiledlayout(1,2, "TileSpacing","compact", "Padding","tight"); 
f.Position = [500 120 800 450];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1 = nexttile; %top left
plot(x_arr_AE(2,1:Tsim), 'Color', col1, "LineWidth", lw2); hold on;
%plot(x_arr_AEperf(2,1:end), 'Color', col2, "LineWidth", lw2, "LineStyle", ls2);
plot(x_arr_Earts(2,1:Tsim), 'Color', col4, "LineWidth", lw4, "LineStyle", ls4);
plot([200; 200], [0.03; 0.13], "LineWidth", lw, "LineStyle","--", 'Color',[0,0,0,grey_alpha]);
ylabel("$$[x]_2$$", "Interpreter", "latex", 'FontSize', fs);
yticks([0.078, 0.082, 0.086, 0.090])
ylim([0.078, 0.092])
grid on;
set(gca,'fontsize',fs-2);
set(gca, "TickLabelInterpreter", "latex");
xlim(ax1, [0, 3*percent*Tsim*1.05])
xlabel("Time", "Interpreter", "latex", 'FontSize', fs);

ax2 = nexttile; % top right
plot(x_arr_AE(2,1:Tsim), 'Color', col1, "LineWidth", lw); hold on;
plot(x_arr_AEperf(2,1:end), 'Color', col2, "LineWidth", lw2, "LineStyle", ls2);
plot(x_arr_Earts(2,1:Tsim), 'Color', col4, "LineWidth", lw4, "LineStyle", ls4);
plot([200; 200], [0.03; 0.13], "LineWidth", lw, "LineStyle","--", 'Color',[0,0,0,grey_alpha]);
xlabel("Time", "Interpreter", "latex", 'FontSize', fs);
%ylabel("$$[x]_2$$", "Interpreter", "latex", 'FontSize', fs);
yticks([0.078, 0.082, 0.086, 0.090])
ylim([0.078, 0.092])
grid on;
set(gca,'fontsize',fs-2);
%ylabel("$$[x]_3$$", "Interpreter", "latex", 'FontSize', fs);
set(gca,'YTickLabel',[]);
set(gca, "TickLabelInterpreter", "latex");
xlim(ax2, [(1-percent)*Tsim*0.997, Tsim*1.003])
leg = legend('AE-MPC', 'E-MPC with $$\theta^\ast$$', 'E-MPC with $$\hat{\theta}$$', 'FontSize', fs, 'Orientation', 'horizontal', 'Interpreter', 'latex');
leg.Layout.Tile = 'north';

linkaxes([ax1 ax2], 'y')
if include_title
    title(t, "Closed-Loop State Evolution", 'FontSize', fs+3, 'Interpreter', 'latex');
end
if saving_fig
    N = input('Horizon length:\n')
    beta = input('Value of beta:\n')
    LMS = 1;
    mu = input('Value of mu:\n')
    comment = input('Additional comments (may be empty):', 's')
    exportgraphics(gcf, ['Results/Figures/States_TI-param' num2str(N) '_' num2str(Tsim) '_beta_' num2str(beta) '_LMS_' num2str(LMS) '_mu_' num2str(mu) '_comment-' comment '.pdf'])
    exportgraphics(gcf, ['Results/Figures/States_TI-param' num2str(N) '_' num2str(Tsim) '_beta_' num2str(beta) '_LMS_' num2str(LMS) '_mu_' num2str(mu) '_comment-' comment '.eps'])
end


%% figure(2): Constant parameter
% Performance comparison, averaged x2 closed-loop evolution and parameter
% estimate's evolution
% close 2
f = figure(2);
lw = 1.75; fs = 16;
f.Position = [0 0 900 900];

plot_type = "averaged_upT"; interv = 10;
% plot_type = "moving_mean";
augmented = 1;
num_av = 5;
val_arr = linspace(0, Tsim, num_av)+1;

run_time = min([length(x_arr_AE(2,:)), length(x_arr_AEperf(2,:)), length(x_arr_Earts(2,:))]);
val_AE = zeros(run_time,1);
val_Earts = zeros(run_time,1);
val_AEperf = zeros(run_time,1);


for jj = 1:run_time
    val_AE(jj) = -sum(x_arr_AE(2,1:jj))/jj;
    val_AEperf(jj) = -sum(x_arr_AEperf(2,1:jj))/jj;
    val_Earts(jj) = -sum(x_arr_Earts(2,1:jj))/jj;
end

run parameter_def.m

if augmented
    t = tiledlayout(4,1, "TileSpacing","compact"); 
end
    ax0 = nexttile;
    plot(x_arr_AE(2,1:Tsim), 'Color', col1, "LineWidth", lw); hold on;
    % plot(x_arr_AEperf(2,1:end), 'Color', col2, "LineWidth", lw2, "LineStyle", ls2);
    plot(x_arr_Earts(2,1:Tsim), 'Color', col4, "LineWidth", lw, "LineStyle", ls4);
    %plot([200; 200], [0.03; 0.13], "LineWidth", lw, "LineStyle","--", 'Color',[0,0,0,0.4]);
    plot([200; 200], [0; 1], "LineWidth", lw+1, "LineStyle","--", 'Color',[0,0,0,grey_alpha]);
    xlabel("Time $$T$$ [-]", "Interpreter", "latex", 'FontSize', fs);
    ylabel("$$[x]_2$$", "Interpreter", "latex", 'FontSize', fs);
    yticks([0.078, 0.082, 0.086, 0.090])
    ylim([0.078, 0.095])
    grid on;
    set(gca,'fontsize',fs-2);
    % set(gca,'XTickLabel',[]);
    set(gca, "TickLabelInterpreter", "latex");
    xlim([0,1000])
    % legend('AE-MPC', 'E-MPC w/ $$\theta^*$$',  'E-MPC w/ $$\hat{\theta}$$', 'FontSize', fs-1, 'Orientation', 'horizontal', 'Interpreter', 'latex', 'Location', 'northeast');
    legend('AE-MPC', 'E-MPC w/ $$\hat{\theta}_0$$', 'FontSize', fs-1, 'Orientation', 'horizontal', 'Interpreter', 'latex', 'Location', 'northeast');

%else
%   t = tiledlayout(3,1, "TileSpacing","compact"); 
%end


% Average of state x2
ax1 = nexttile; hold on; grid on;
av_length = 50; av_opt = "shrink";

switch plot_type
    case "moving_mean"
        plot(-movmean(x_arr_AE(2,:),av_length, "Endpoints", av_opt), 'Color', col1, "LineWidth", lw, 'LineStyle', '-');
        plot(-movmean(x_arr_Earts(2,:),av_length, "Endpoints", av_opt), 'Color', col4, "LineWidth", lw, 'LineStyle', ls4);
         plot(-movmean(x_arr_AEperf(2,:),av_length, "Endpoints", av_opt), 'Color', col2, "LineWidth", lw, 'LineStyle', ls2);

    case "averaged_upT"
        plot(1:interv:run_time, val_AE(1:interv:end), 'Color', col1, "LineWidth", lw, 'LineStyle', '-');
        plot(1:interv:run_time, val_Earts(1:interv:end), 'Color', col4, "LineWidth", lw, 'LineStyle', ls4);
        plot(1:interv:run_time, val_AEperf(1:interv:end), 'Color', col2, "LineWidth", lw, 'LineStyle', ls2);
    otherwise
        error("Wrong case.")
end

plot([0,Tsim], ones(2,1)*-0.084649344349529, "LineWidth", lw, "LineStyle","-.", 'Color',[0,0,0,grey_alpha]), 
plot([200; 200], [-0.1; 0], "LineWidth", lw+1, "LineStyle","--", 'Color',[0,0,0,grey_alpha]);
ylim([-0.0875, -0.0845])
set(gca,'Box','on');
set(gca,'fontsize',fs-2);
% set(gca,'XTickLabel',[]);
xlabel("Time $$T$$ [-]", "Interpreter", "latex", 'FontSize', fs);
set(gca, "TickLabelInterpreter", "latex");
ylabel("$$\sum_{k=0}^T\ell(x_k,u_k)/T$$", "Interpreter", "latex", 'FontSize', fs-1);
legend('AE-MPC', 'E-MPC w/ $$\hat{\theta}_0$$', 'E-MPC w/ $$\theta_k$$',  'FontSize', fs-1, 'Orientation', 'horizontal', 'Interpreter', 'latex', 'Location', 'southeast');

% Parameter estimate 1
ax2 = nexttile; hold on; grid on;
plot(theta_arr_AE(1,1:end), "LineWidth", lw);
%theta_true lines
true_p11 = plot([0;Tsim], theta_true(1)*[1;1], 'Color', coltrue, 'LineStyle', lstrue, "LineWidth", lw);
%Theta set
plot([0;Tsim], max(Theta_0.V(:,1))*[1;1], 'Color', colcon, 'LineStyle', lscon,  "LineWidth", lw); 
plot([0;Tsim], min(Theta_0.V(:,1))*[1;1], 'Color', colcon, 'LineStyle', lscon,  "LineWidth", lw); 
plot([200; 200], [0; 1.5], "LineWidth", lw+1, "LineStyle","--", 'Color',[0,0,0,grey_alpha]);
ylim([0.995*min(Theta_0.V(:,1)), 1.005*max(Theta_0.V(:,1))])
set(gca,'Box','on');
set(gca,'fontsize',fs-2);
set(gca,'XTickLabel',[]);
set(gca, "TickLabelInterpreter", "latex");
legend('$$[\hat{\theta}]_1$$', '$$[\theta]_1$$', '$$\Theta$$' , 'Interpreter', 'latex', 'FontSize', fs, "Orientation", "horizontal", "Location", "southeast"); 
% xlabel("Time [-]", "Interpreter", "latex", 'FontSize', fs);
ylabel("$$[\theta]_1$$", "Interpreter", "latex",  'FontSize', fs); 


% Parameter estimate 2
ax3 = nexttile; hold on; grid on;
plot(theta_arr_AE(2,1:end), "LineWidth", lw); 
%theta_true lines
true_p21 = plot([0;Tsim], theta_true(2)*[1;1], 'Color', coltrue, 'LineStyle', lstrue, "LineWidth", lw);
%Theta set
plot([0;Tsim], max(Theta_0.V(:,2))*[1;1], 'Color', colcon,  'LineStyle', lscon, "LineWidth", lw);
plot([0;Tsim], min(Theta_0.V(:,2))*[1;1], 'Color', colcon,  'LineStyle', lscon,  "LineWidth", lw);
plot([200; 200], [0.5; 1.5], "LineWidth", lw+1, "LineStyle","--", 'Color',[0,0,0,grey_alpha]);
xlabel("Time $$T$$ [-]", "Interpreter", "latex", 'FontSize', fs); 
ylabel("$$[\theta]_2$$", "Interpreter", "latex",  'FontSize', fs); 
ylim([0.995*min(Theta_0.V(:,2)), 1.005*max(Theta_0.V(:,2))])
set(gca,'Box','on');
set(gca,'fontsize',fs-2);
set(gca, "TickLabelInterpreter", "latex");
legend('$$[\hat{\theta}]_2$$', '$$[\theta]_2$$', '$$\Theta$$' , 'Interpreter', 'latex', 'FontSize', fs, "Orientation", "horizontal"); 


% if augmented
    linkaxes([ax1 ax2 ax3], 'x')
%else
%    linkaxes([ax1 ax2 ax3], 'x')
% end
xlim([0, Tsim*1.01])
if include_title
    title(t, "Performance of AE-MPC compared with E-MPC", 'FontSize', fs+2, 'Interpreter', 'latex')
end
if saving_fig
    disp("Please input data for naming the figure as follows:")
    N = input('Horizon length:\n')
    beta = input('Value of beta:\n')
    LMS = 1;
    mu = input('Value of mu:\n')
    comment = input('Additional comments (may be empty):', 's')
    exportgraphics(gcf, ['Results/Figures/Perf_TIparam--' num2str(N) '_' num2str(Tsim) '_beta_' num2str(beta) '_LMS_' num2str(LMS) '_mu_' num2str(mu) '_comment-' comment '.pdf'])
    % exportgraphics(gcf, ['Results/Figures/Perf_TIparam--' num2str(N) '_' num2str(Tsim) '_beta_' num2str(beta) '_LMS_' num2str(LMS) '_mu_' num2str(mu) '_comment-' comment '.eps'])
end


%% figure(6): Time-varying parameter
% Performance comparison, averaged x2 closed-loop evolution and parameter
% estimate's evolution
TV = true;
run parameter_def.m

%%% load the correct data
load("Results/Simulations/AE_scaled--N_25-Tsim_3000-mu_15-h_0.025-LMS_1-TV_1.mat", "x_arr", "theta_arr", "Tsim")
x_arr_AE = x_arr;
theta_arr_AE_TVp = theta_arr;

% load("Results/055/LMI_designed/TV/E_wrongTV--N_25-Tsim_1000-mu_150-h_0.025-LMS_0-TV_1.mat", "x_arr", "theta_arr")
load("Results/Simulations/E_true--N_25-Tsim_3000-mu_15-h_0.025-LMS_0-TV_1.mat", "x_arr", "theta_arr")
x_arr_AEperf = x_arr;
theta_arr_AEperf = theta_arr;

% load("Results/055/LMI_designed/TV/E_trueTV2--N_25-Tsim_4000-mu_150-h_0.025-LMS_0-TV_1.mat", "x_arr", "theta_arr")
load("Results/Simulations/E_wrong--N_25-Tsim_3000-mu_15-h_0.025-LMS_0-TV_1.mat", "x_arr", "theta_arr")
x_arr_Earts = x_arr;
theta_arr_Earts = theta_arr;

% preliminary things
grey_alpha = 0.5; lw = 1.75; fs = 16;
run_time = min([length(x_arr_AE(2,:)), length(x_arr_AEperf(2,:))]);%, length(x_arr_Earts(2,:))]);
% rampval = [0.0000004;-0.000007];
rampval = [0.00000025;-0.000004]; 
theta_true = theta_true.*[0.995;1.01]; % alters the true parameter value for the time-varying parameter case, must match the value in AE_MPC_journal.m!
sw_time = [0,Tsim];
ramp_end = [rampval(1)*(sw_time(2)-sw_time(1)); rampval(2)*(sw_time(2)-sw_time(1))];



figure(6); close;
f = figure(6);
f.Position = [300 60 900 725];
t = tiledlayout(3,1, "TileSpacing","compact", "Padding","loose"); 

num_av = 5;
val_arr = linspace(0, Tsim, num_av)+1;
val_AE = zeros(num_av-1,1);
val_Earts = zeros(num_av-1,1);
val_AEperf = zeros(num_av-1,1);
% sw_time = [1500; 2500]; 
sw_time = [0; Tsim]; 

% Set the default font 
% set(groot, 'defaultAxesFontName', 'Arial'); 
% set(groot, 'defaultTextFontName', 'Arial');

for jj = 1:run_time
    val_AE(jj) = -sum(x_arr_AE(2,1:jj))/jj;
    val_AEperf(jj) = -sum(x_arr_AEperf(2,1:jj))/jj;
    val_Earts(jj) = -sum(x_arr_Earts(2,1:jj))/jj;
end

% Average of state x2
ax1 = nexttile; hold on; grid on;
av_length = 150; av_opt = "shrink";
plot(val_AE, 'Color', col1, "LineWidth", lw*1.5, 'LineStyle', '-');
plot(val_AEperf, 'Color', col2, "LineWidth", lw*1.5, 'LineStyle', ':');
plot(val_Earts, 'Color', col4, "LineWidth", lw*1.5, 'LineStyle', '--');
plot([200; 200], [-1; 0], "LineWidth", lw+1, "LineStyle","--", 'Color',[0,0,0,grey_alpha], "HandleVisibility","off");
% xlabel("Time [-]", "Interpreter", "latex", 'FontSize', fs); 
ylabel("$$\sum_{k=0}^T\ell(x_k,u_k)/T$$", "Interpreter", "latex",  'FontSize', fs); 
ylim([-0.0862, -0.0838]);
set(gca,"YTick", [-0.086,-0.085,-0.084]);
set(gca,'Box','on');
set(gca,'fontsize',fs-2);
set(gca, "TickLabelInterpreter", "latex");
set(gca,'XTickLabel',[]);
legend('AE-MPC', 'E-MPC w/ $$\theta_0$$',  'E-MPC w/ $$\hat{\theta}_k$$', 'FontSize', fs-1, 'Orientation', 'horizontal', 'Interpreter', 'latex', "Location", "southeast");
%set(gca, 'FontName','Times New Roman');

% Parameter estimate 1
ax2 = nexttile;
plot(theta_arr_AE_TVp(1,1:end), "LineWidth", lw); 
hold on; grid on;
%theta_true lines
plot([0;sw_time(1)], theta_true(1)*[1;1], 'Color', coltrue, 'LineStyle', lstrue, "LineWidth", lw);
plot([sw_time], [theta_true(1);theta_true(1)+ramp_end(1)], 'Color', coltrue, 'LineStyle', lstrue, "LineWidth", lw, "HandleVisibility","off");
plot([sw_time(2); Tsim], (theta_true(1)+ramp_end(1))*[1;1], 'Color', coltrue, 'LineStyle', lstrue, "LineWidth", lw, "HandleVisibility","off");
plot([200; 200], [0; 1.5], "LineWidth", lw+1, "LineStyle","--", 'Color',[0,0,0,grey_alpha], "HandleVisibility","off");
%Theta set
plot([0;Tsim], max(Theta_0.V(:,1))*[1;1], 'Color', colcon, 'LineStyle', lscon,  "LineWidth", lw); 
plot([0;Tsim], min(Theta_0.V(:,1))*[1;1], 'Color', colcon, 'LineStyle', lscon,  "LineWidth", lw); 
ylim([0.995*min(Theta_0.V(:,1)), 1.005*max(Theta_0.V(:,1))])
legend('$$[\hat{\theta}_k]_1$$', '$$[\theta_k]_1$$', '$$\Theta$$' , 'Interpreter', 'latex', 'FontSize', fs, "Orientation", "horizontal", "Location", "southeast"); 
% xlabel("Time [-]", "Interpreter", "latex", 'FontSize', fs) 
ylabel("$$[\theta]_1$$", "Interpreter", "latex",  'FontSize', fs); 
set(gca,'Box','on');
set(gca,'fontsize',fs-2);
set(gca,'XTickLabel',[]);
set(gca, "TickLabelInterpreter", "latex");
% set(gca, 'FontName','Arial');


% Parameter estimate 2
ax3 = nexttile;
plot(theta_arr_AE_TVp(2,1:end), "LineWidth", lw); 
hold on; grid on;
%theta_true lines
plot([0;sw_time(1)], theta_true(2)*[1;1], 'Color', coltrue, 'LineStyle', lstrue, "LineWidth", lw);
plot([sw_time], [theta_true(2);theta_true(2)+ramp_end(2)], 'Color', coltrue, 'LineStyle', lstrue, "LineWidth", lw, "HandleVisibility","off");
plot([sw_time(2); Tsim], (theta_true(2)+ramp_end(2))*[1;1], 'Color', coltrue, 'LineStyle', lstrue, "LineWidth", lw, "HandleVisibility","off");
plot([200; 200], [0; 1.5], "LineWidth", lw+1, "LineStyle","--", 'Color',[0,0,0,grey_alpha], "HandleVisibility", "off");
%Theta set
plot([0;Tsim], max(Theta_0.V(:,2))*[1;1], 'Color', colcon,  'LineStyle', lscon, "LineWidth", lw);
plot([0;Tsim], min(Theta_0.V(:,2))*[1;1], 'Color', colcon,  'LineStyle', lscon,  "LineWidth", lw);
xlabel("Time T [-]", "Interpreter", "latex", 'FontSize', fs); 
ylabel("$$[\theta]_2$$", "Interpreter", "latex",  'FontSize', fs); 
ylim([0.995*min(Theta_0.V(:,2)), 1.005*max(Theta_0.V(:,2))])
leg = legend("$$[\hat{\theta}_k]_2$$", "$$[\theta_k]_2$$", "$$\Theta$$", 'Interpreter', 'latex', 'FontSize', fs, "Orientation", "horizontal"); 
set(gca,'Box','on');
set(gca,'fontsize',fs-2);
set(gca, "TickLabelInterpreter", "latex");
% set(gca, 'FontName','Times New Roman');


linkaxes([ax1 ax2 ax3], 'x')
xlim([0, Tsim*1.01])
if include_title
    title(t, "Performance Comparison for Time-Varying Parameters", 'FontSize', fs+2, 'Interpreter', 'latex');
end

%% figure(63): State x2 at beginning and end of the simulation
f = figure(13); t = tiledlayout(1,2, "TileSpacing","compact", "Padding","tight"); 
f.Position = [500 120 800 450];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1 = nexttile; %top left
plot(x_arr_AE(2,1:Tsim), 'Color', col1, "LineWidth", lw2); hold on;
%plot(x_arr_AEperf(2,1:end), 'Color', col2, "LineWidth", lw2, "LineStyle", ls2);
plot(x_arr_Earts(2,1:Tsim), 'Color', col4, "LineWidth", lw4, "LineStyle", ls4);
plot([200; 200], [0.03; 0.13], "LineWidth", lw, "LineStyle","--", 'Color',[0,0,0,grey_alpha]);
ylabel("$$[x]_2$$", "Interpreter", "latex", 'FontSize', fs);
yticks([0.078, 0.082, 0.086, 0.090])
ylim([0.078, 0.092])
grid on;
set(gca,'fontsize',fs-2);
set(gca, "TickLabelInterpreter", "latex");
xlim(ax1, [0, 3*percent*Tsim*1.05])
xlabel("Time", "Interpreter", "latex", 'FontSize', fs);

ax2 = nexttile; % top right
plot(x_arr_AE(2,1:Tsim), 'Color', col1, "LineWidth", lw); hold on;
plot(x_arr_AEperf(2,1:end), 'Color', col2, "LineWidth", lw2, "LineStyle", ls2);
plot(x_arr_Earts(2,1:Tsim), 'Color', col4, "LineWidth", lw4, "LineStyle", ls4);
plot([200; 200], [0.03; 0.13], "LineWidth", lw, "LineStyle","--", 'Color',[0,0,0,grey_alpha]);
xlabel("Time", "Interpreter", "latex", 'FontSize', fs);
%ylabel("$$[x]_2$$", "Interpreter", "latex", 'FontSize', fs);
yticks([0.078, 0.082, 0.086, 0.090])
ylim([0.078, 0.092])
grid on;
set(gca,'fontsize',fs-2);
%ylabel("$$[x]_3$$", "Interpreter", "latex", 'FontSize', fs);
set(gca,'YTickLabel',[]);
set(gca, "TickLabelInterpreter", "latex");
xlim(ax2, [(1-percent)*Tsim*0.997, Tsim*1.003])
leg = legend('AE-MPC', 'E-MPC with $$\theta^\ast$$', 'E-MPC with $$\hat{\theta}$$', 'FontSize', fs, 'Orientation', 'horizontal', 'Interpreter', 'latex');
leg.Layout.Tile = 'north';

linkaxes([ax1 ax2], 'y')
if include_title
    title(t, "Closed-Loop State Evolution", 'FontSize', fs+3, 'Interpreter', 'latex');
end
if saving_fig
    N = input('Horizon length:\n')
    beta = input('Value of beta:\n')
    LMS = 1;
    mu = input('Value of mu:\n')
    comment = input('Additional comments (may be empty):', 's')
    exportgraphics(gcf, ['Results/Figures/States_TI-param' num2str(N) '_' num2str(Tsim) '_beta_' num2str(beta) '_LMS_' num2str(LMS) '_mu_' num2str(mu) '_comment-' comment '.pdf'])
    exportgraphics(gcf, ['Results/Figures/States_TI-param' num2str(N) '_' num2str(Tsim) '_beta_' num2str(beta) '_LMS_' num2str(LMS) '_mu_' num2str(mu) '_comment-' comment '.eps'])
end



if saving_fig
    disp("Please input data for naming the figure as follows:")
    N = input('Horizon length:\n')
    beta = input('Value of beta:\n')
    LMS = 1;
    mu = input('Value of mu:\n')
    comment = input('Additional comments (may be empty):', 's')
    figure(6);
    exportgraphics(gcf, ['Results/Figures/Perf_TVparam--' num2str(N) '_' num2str(Tsim) '_beta_' num2str(beta) '_LMS_' num2str(LMS) '_mu_' num2str(mu) '_comment-' comment '.pdf'])
    % exportgraphics(gcf, ['Results/Figures/Perf_TIparam--' num2str(N) '_' num2str(Tsim) '_beta_' num2str(beta) '_LMS_' num2str(LMS) '_mu_' num2str(mu) '_comment-' comment '.eps'])
end

