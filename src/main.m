%% main.m
% This script starts the offline computations and then the simulation.
clear all; clc; close all;

diary_active = true; % create log file with console output if set to TRUE, no log file if set to FALSE

if diary_active
    diary_file = string(datetime, "yyyy-MM-dd_hh-mm-ss");
    diary(diary_file);
    diary on;
    disp("Starting main.m at time "+ string(datetime, "yyyy-MM-dd_hh:mm:ss"));
end

%% Offline computations
% dynamics, parameters and noise are defined in syst.m and parameter_def.m
run offline_pipeline.m

%% Simulations
% By default, the just-in-time compilation of CasADi's solver is disabled.
% Just-in-time compilation can be used to reduce the online computation time.
% To enable this feature, set the JIT boolean to true. For Windows users: Be sure to open Matlab from the 
% command prompt of the "x64 Native Tools Command Propmt for VS" (usually by typing "matlab" + ENTER ). 
% For Linux/Mac users: Be sure to have "gcc" installed. For more details see casadi.org and 
% https://github.com/casadi/casadi/wiki/FAQ:-how-to-perform-jit-for-function-evaluations-of-my-optimization-problem%3F
JIT = false;

%%%%%%
% time-invariant parameters
N = 25;         % prediction horizon
beta = 100;     % weight of artifical reference cost
mu = 15;        % parameter update gain (least-mean squares estimator)
Tsim = 2500;      % total simulation time

AE_MPC("AEMPC", N, Tsim, beta, mu, false, JIT);       % Simulation of AE-MPC
AE_MPC("EMPC_estimate", N, Tsim, beta, mu, false, JIT);       % Simulation of E-MPC (no adaptation)
AE_MPC("EMPC_perfect", N, Tsim, beta, mu, false, JIT);       % Simulation of E-MPC (perfect model knowledge)
disp("END of time-invariant simulations.")

%%%%%%
% time-varying parameters
N = 25;         % prediction horizon
beta = 100;     % weight of artifical reference cost
mu = 15;        % parameter update gain (least-mean squares estimator)
Tsim = 3000;    % total simulation time

AE_MPC("AEMPC", N, Tsim, beta, mu, true, JIT);       % Simulation of AE-MPC
AE_MPC("EMPC_estimate", N, Tsim, beta, mu, true, JIT);       % Simulation of E-MPC (no adaptation)
AE_MPC("EMPC_perfect", N, Tsim, beta, mu, true, JIT);       % Simulation of E-MPC (perfect model knowledge)
disp("END of timevarying simulations.")

%% Plotting
% For saving figures, please follow the console instructions.
warning("If you changed the simulation parameters (horizon, Tsim, beta, mu), you must change the filenames in the plotting script accordingly.")
run plotting_figures.m

if diary_active
    disp("END of main.m at time " + string(datetime, "yyyy-MM-dd_hh:mm:ss"));
    diary off;
else
    disp("END of main.m")
end

