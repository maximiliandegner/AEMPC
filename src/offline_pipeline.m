%% offline_pipeline.m
% Option 1: fix values for w_bar, rho, L_w and solve one SDP
% Option 2: solve the SDP for a range of rho and L_w values

%% clear everything and get model parameters and uncertainty
% clear all; close all; clc;
t = tic;

% adjust model parameters and uncertainty if desired --> change the values in the file below
addpath("Offline\");
theta_true = syst.theta_true;
theta_nom = syst.theta_nom;
run parameter_def.m
a_1 = syst.scl(1); a_2 = syst.scl(2);      % scaling factors for parameter in model
option = 1

if option == 1
%% Option 1: fix w_bar, rho, L_w and solve SDP
w_max = 0.002;      % max. uncerainty over \Theta and \mathbb{D}
w_bar = w_max;
rho = 0.97664;      % contraction rate
Lw = 0.0074444;     % Lipschitz constant
h = 0.025;          % discretization time
n_steps = 3;

run constraint_def.m
run opt_steady_offline.m
load('opt_steady.mat', "r_set"); equilib = full(r_set);


% preliminary definitions
import casadi.*  
nonconst = 0;

% combined state and input constraints + evaluation at desired
% steady-state
H_xu = [eye(n+m); -eye(n+m)];
h_xu = [b_x(1:n); b_u(1:m); b_x(n+1:end); b_u(m+1:end)];
g_val = H_xu*equilib - h_xu;

% nonlinear, parameter-dependent part of the dynamics
xx = MX.sym('xx', 3,1);
delta = 0.55;
G = h* [-xx(1)^2*exp(-1/xx(3)), -xx(1)*exp(-delta/xx(3));
    xx(1)^2*exp(-1/xx(3))  ,  MX(0);
    MX(0)                  ,  MX(0)];
G_f = Function('G_f', {xx}, {G});
G_1x = Function('G_1x', {xx}, {jacobian(G(:,1), xx)});
G_2x = Function('G_2x', {xx}, {jacobian(G(:,2), xx)});

% creates the LMIs and solves the SDP
addpath("Offline\");
run solve_offline_SDP.m
        % creates the variable 'solved' to indicate solve status

% obtain solution if solved successfully
if solved
    disp("Solution found.")
    Y=value(Y);  Y_0=value(Y_0);
    X=value(X);  X_0=value(X_0);
    P = value(inv(X));
    K = value(Y*P);
else
    disp("Not solved.")
end

% save values to "Offline/grid_disc.mat"
save("Offline/grid_disc.mat", 'Y', 'Y_0', 'X_0','X','h', "rho", "K", "P", "w_max", 'Lw');
disp("End of Option 1.")

elseif option == 2
%% Option 2: or fix w_max, specify range of L_w and rho and solve the SDP repeatedly
h = 0.025; 
n_steps = 3;
RHO = linspace(0.97, 0.9999, 10);
LW_range = linspace(0.001, 0.03, 10);
w_max = 0.002; 

% solve the SDP for a range of values
% interaction is required to select the solution
[P, K, w_bar, rho, c_j, Lw] = LMI_improved(h, n_steps, RHO, LW_range, w_max);

disp("End of Option 2.")


else
    error("Select an option (1 or 2).")
end
%% compute terminal cost
run compute_terminal_cost.m

toc(t)

%% display results in console
disp("Tube results:")
P
K
w_bar
rho
Lw


disp("Terminal cost:")
H
alpha = M
P_f

%% next: Simulation
% --> AE_MPC_journal(N, Tsim, beta, mu, LMS, TV)
%
% e.g.:
% N = 25;       % prediction horizon
% Tsim = 100;   % simulation time
% beta = 100;   % weight on artificial reference cost term
% mu = 150;     % LMS update gain
% LMS = 1;      % LMS on: 1, off: 0
% TV = 0;       % time-varying parameter: 1; constant parameter: 0

% See script AE_MPC_journal for more info
