%% main.m
% This script starts the offline computations and then the simulation.
clear all; clc; close all;

%% Offline computations
run offline_pipeline.m

%% Simulations
% By default, the just-in-time compilation of CasADi's solver is disabled.
% To enable this feature, set the JIT boolean to true. For more details see
% casadi.org and https://github.com/casadi/casadi/wiki/FAQ:-how-to-perform-jit-for-function-evaluations-of-my-optimization-problem%3F
JIT = false;

%%%%%%
% time-invariant parameters
N = 25;         % prediction horizon
beta = 100;     % weight of artifical reference cost
mu = 15;        % parameter update gain (least-mean squares estimator)
Tsim = 10;      % total simulation time

AE_MPC_journal("AEMPC", N, Tsim, beta, mu, false, JIT);       % Simulation of AE-MPC
AE_MPC_journal("EMPC_estimate", N, Tsim, beta, mu, false, JIT);       % Simulation of E-MPC (no adaptation)
AE_MPC_journal("EMPC_perfect", N, Tsim, beta, mu, false, JIT);       % Simulation of E-MPC (perfect model knowledge)
disp("END of time-invariant simulations.")

%%%%%%
% time-varying parameters
N = 25;         % prediction horizon
beta = 100;     % weight of artifical reference cost
mu = 15;        % parameter update gain (least-mean squares estimator)
Tsim = 3000;    % total simulation time

AE_MPC_journal("AEMPC", N, Tsim, beta, mu, true, JIT);       % Simulation of AE-MPC
AE_MPC_journal("EMPC_estimate", N, Tsim, beta, mu, true, JIT);       % Simulation of E-MPC (no adaptation)
AE_MPC_journal("EMPC_perfect", N, Tsim, beta, mu, true, JIT);       % Simulation of E-MPC (perfect model knowledge)
disp("END of timevarying simulations.")

%% Plotting
% For saving figures, please follow the console instructions.
warning("If you changed the simulation parameters (horizon, Tsim, beta, mu), you must change the filenames in the plotting script accordingly.")
run figures_journal.m
disp("END of main.m")
