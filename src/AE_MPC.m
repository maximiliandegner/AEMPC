% clear all
function AE_MPC(type, N, Tsim, beta, mu, TV, JIT_boolean)
 clc
% type      String, one of the following options: "AEMPC", "EMPC_perfect", or "EMPC_estimate"
% N         Prediction horizon
% Tsim      Simulation length
% beta      Weight of artifical reference cost term (equilibrium cost)
% mu        LMS parameter update gain
% JIT_boolean   Switch just-in-time compilation (CasADi) of Solver on
                % (true) or off (false); default: off


%% Input parsing
switch nargin
    case 0
        type = "AEMPC";
        N = 15;
        Tsim = 200;
        beta = 100;
        mu = 15;
        TV = false;
        JIT_boolean = false;
    case 6
        JIT_boolean = false;
    case 7

    otherwise
        error("Number of input arguments incorrect.")
end
if ~(type == "AEMPC" || type == "EMPC_estimate" || type == "EMPC_perfect")
    error("Please enter an allowed AE-MPC or E-MPC type.")
end

if type == "AEMPC"
    LMS = true; warning("Switched LMS <ON> due to AE-MPC mode.")
elseif  contains(type,"EMPC")
    LMS = false; warning("Switched LMS <OFF> due to E-MPC mode.")
end

if JIT_boolean
    warning("Please ensure that MATLAB has access to the defined compiler, see solver options.")
else
    warning("Running <without> just-in-time compilation.")
end


%% import and system parameters
import casadi.* 
n = syst.n; % state dimension
m = syst.m; % input dimension 
q = syst.q; % parameter dimension
theta_true = syst.theta_true;
theta_nom = syst.theta_nom;

% load parameter defintion --> gives Theta_0, W
run parameter_def.m;

if TV
    % alters the true parameter value for the time-varying parameter case, must match the value in figures_journal.m / Figure(6)!
    theta_true = theta_true.*[0.995;1.01];  
end
load("Offline/grid_disc.mat", "K", "P", "rho", "h", "w_max", "Lw")
fprintf("Input arguments: N= %d, Tsim= %d, beta= %d, mu= %d, h= %d\n",[N, Tsim, beta, mu, h]);
pause(2);

w_bar = w_max;

% display feedback and tube parameters in console
K
P
rho
s_bar = w_bar/(1-rho-Lw)
load("opt_steady.mat", "r_set")     % loads optimal steady-state for the true parameter

rng(2);
con_obj.uncbound_nlpinit = true;

%% Struct setup: Constraints, 

% simulation time start
t_0 = 0.0;

%constraints (x_min, etc. stem from parameter_def.m)
con_obj.ub_x = [x_max;x_max;x_max];
con_obj.lb_x = [x_min;x_min;x_min];
con_obj.ub_u = u_max;
con_obj.lb_u = u_min;
con_obj.n = n; con_obj.m = m;
con_obj.lambda = inf;
con_obj.solver_lambda = [];
con_obj.solver_uncertbound = [];
con_obj.Theta = Theta_0;
con_obj.theta_nom = theta_nom;
con_obj.W = W;
con_obj.P = P;

% referenc constraints
xr_min = x_min; xr_max = x_max;
ur_min=0.049;
ur_max=0.449;

A_x = [eye(n); -eye(n)];    b_x = [xr_max;xr_max;xr_max;  -xr_min;-xr_min;-xr_min];
A_u = [eye(m); -eye(m)];    b_u = [ur_max; -ur_min];


%% Optimization variable definition
opti = casadi.Opti();

x = opti.variable((N+1)*n,1);
u = opti.variable(N*m,1);
xr = opti.variable(n,1);
ur = opti.variable(m,1);
zr = opti.variable(n,1);
vr = opti.variable(m,1);
z = opti.variable((N+1)*n,1);
v = opti.variable(N*m,1);
s = opti.variable(N+1,1);
w = opti.variable(N,1);
p = opti.variable(n,1);

x_k             = opti.parameter(n, 1);
theta_LMS       = opti.parameter(q, 1);
lambda_old      = opti.parameter(q, 1);

%% Set initial parameter estimate and define objective function
switch type
    case "AEMPC"
        theta_est = [0.09803; 1.02];
        disp("INFO:   Running AE-MPC mode.")
    case "EMPC_estimate"
        theta_est = [0.09803; 1.02];
        disp(['INFO:   Starting E-MPC mode with estimate = ' num2str(theta_est(1)) ', ' num2str(theta_est(2)) '.'])
    case "EMPC_perfect"
        theta_est = theta_true;
        disp(['INFO:   Running E-MPC mode with true = ' num2str(theta_est(1)) ', ' num2str(theta_est(2)) '.'])
end
con_obj.theta_nom = theta_nom;
con_obj.theta_true = theta_true;

% loading terminal cost
load("Offline\terminal_cost_P_f.mat", "P_f")
run compute_constraint_tightening.m;

% definition of the objective
objective = sum(-x(2:n:end-n)) + p'*(x(end-n+1:end)-xr) + (x(end-n+1:end)-xr)'*P_f*(x(end-n+1:end)-xr) + beta*(-xr(2));


%% Define constraints
concheck = []; % array used for debugging and finding violated constraints
opti.minimize(  objective   );
for ii = 1:N
    % state constraints, input constraints
    opti.subject_to( A_x*z((ii-1)*n+1:ii*n) - b_x + c(1:size(A_x,1))*s(ii) <= zeros(size(A_x,1),1) ); concheck = [concheck; -ones(6,1)];
    opti.subject_to( A_u*v((ii-1)*m+1:ii*m) - b_u + c(size(A_x,1)+1:end)*s(ii) <= zeros(size(A_u,1),1) );concheck = [concheck; -ones(2,1)];
    
    % dynamics (nominal and LMS)
    opti.subject_to( u((ii-1)*m+1:ii*m) -( K*(x((ii-1)*n+1:ii*n) - z((ii-1)*n+1:ii*n)) + v((ii-1)*m+1:ii*m)  ) == zeros(m,1) );concheck = [concheck; zeros(m,1)];
    opti.subject_to( z((ii)*n+1:(ii+1)*n) -( syst.dynamic(z((ii-1)*n+1:ii*n), v((ii-1)*m+1:ii*m), h, theta_nom)) == zeros(n,1) );concheck = [concheck; zeros(n,1)];
    opti.subject_to( x((ii)*n+1:(ii+1)*n) -( syst.dynamic(x((ii-1)*n+1:ii*n), u((ii-1)*m+1:ii*m), h, theta_LMS)) == zeros(n,1) );concheck = [concheck; zeros(n,1)];
    
    % tube dynamics
    opti.subject_to( s(ii+1) - ( s(ii)*rho + w(ii)) == 0) ; concheck = [concheck; 0];
    

    temp_w = uncertainty_fcn(z((ii-1)*n+1:ii*n), v((ii-1)*m+1:ii*m), h, theta_nom, Theta_0, W, P);   
    max_temp = MX(0);
    for iii = 1:size(temp_w)
        max_temp = fmax(max_temp, temp_w(iii));
    end
    opti.subject_to( w(ii) -max_temp-Lw*s(ii) >= 0); concheck = [concheck; 1];

end
opti.subject_to( -s <= zeros(N+1,1) ); concheck = [concheck; -ones(N+1,1)];

% terminal state constraints for trajectories
ii = N + 1;
opti.subject_to( z((ii-2)*n+1:(ii-1)*n) - z((ii-1)*n+1:ii*n) == zeros(n,1)); concheck = [concheck; zeros(n,1)];

% steady states
opti.subject_to( xr - syst.dynamic(xr, ur, h, theta_LMS) == zeros(n,1));    concheck = [concheck; zeros(n,1)];  %full(r_set(1:3)) ); %
opti.subject_to( zr - z(end-2:end) == zeros(n,1));          concheck = [concheck; zeros(n,1)];
opti.subject_to( ur -( K*(xr-zr) + v(end)) == zeros(m,1) ); concheck = [concheck; zeros(m,1)];
    % terminal state constraints for steady states
opti.subject_to( A_x*z(end-2:end) - b_x + c(1:size(A_x,1))*s(ii) <= zeros(size(A_x,1),1)  ); concheck = [concheck; -ones(6,1)];

% anti performance-degradation constraint
obj.solver_lambda = [];
u_s = full(r_set(4)); x_s = full(r_set(1:3));
opti.set_value(lambda_old,  opt_steady_online(n,m, theta_est,u_s , x_s, true, obj)  );
opti.subject_to( -xr(2) <= lambda_old+1e-3 ); concheck = [concheck; -1];

% initial constraint
opti.subject_to( x(1:n) - x_k == zeros(n,1) );                      concheck = [concheck; zeros(n,1)];
opti.subject_to( (x(1:n)-z(1:n))'*P*(x(1:n)-z(1:n)) <= s(1)^2 );    concheck = [concheck; -1];

% terminal cost linear part
[A,B] = syst.getA_d(z(end-2:end), v(end), h, theta_LMS);
opti.subject_to(  (eye(n) - (A+B*K))' * p - [0;-1;0] == zeros(n,1)); concheck = [concheck; zeros(n,1)];

% terminal constraint, tube
opti.subject_to( s(end) - s(end-1) <= 0); concheck = [concheck; -1];


%% Initializing storage arrays
x_arr = zeros(n,Tsim);              x_arr(1:n,1) = full(r_set(1:n));
z_arr = zeros(n,Tsim);              z_arr(1:n,1:2) = repmat(full(r_set(1:n)),1,2);
xr_arr = zeros(n,Tsim);             xr_arr(1:n,1:2) = repmat(full(r_set(1:n)),1,2);
zr_arr = zeros(n,Tsim);             zr_arr(1:n,1) = full(r_set(1:n));
u_arr = zeros(m,Tsim);              u_arr(1:m,1) = full(r_set(n+1:n+m));
v_arr = zeros(m,Tsim);              v_arr(1:m,1) = full(r_set(n+1:n+m));
theta_arr = zeros(2,Tsim);          theta_arr(:,1:2) = repmat(theta_est, 1,2);
w_arr = zeros(n,Tsim);
Tsolve_arr = zeros(1,Tsim);

x_sol = x_arr(1:n,1);
u_sol = full(r_set(n+1:n+m));
x_next = syst.dynamic(x_sol, u_sol, h, theta_true);
x_arr(1:n, 2) = x_next;


%% Solver options and initialize solver
opts.ipopt.print_level = 0;
% opts.ipopt.derivative_test = 'second-order';
opts.print_time = 0;
    % JUST IN TIME Compilation
opts.jit = JIT_boolean;
opts.compiler = 'shell';
% opts.jit_options.flags = {'-O3'}; % does not work with VS C++ Build Tools
% opts.jit_options.flags = {'/O2'}; % compiler options; '/O1' for small
                                    % code, '/O2' for fast code (/O2 elongates compile time significantly...)
opts.jit_options.verbose = true;
    % SOLVER Initialization
opti.solver('ipopt', opts);

%% setup LMS estimator struct
LMS_params.const.mu             = mu;
LMS_params.const.theta_true     = theta_true;
LMS_params.sys.n                = n;
LMS_params.const.theta_nom      = theta_nom;

%% ONLINE LOOP
for t = 2:Tsim

    %%% SWITCH LMS ON AFTER 200 STEPS
    if LMS && t>= 200
        theta_est = param_estim(x_arr(1:n, t-1), u_arr(1:m, t-1), x_arr(:,t), h, theta_est, Theta_0, LMS_params)
    elseif ~LMS && type == "EMPC_estimate" 
        theta_est =  theta_est; 
    elseif ~LMS && type == "EMPC_perfect" 
        theta_est =  theta_true; 
    end

    opti.set_value(theta_LMS, theta_est);
    opti.set_value(x_k, x_next);

    if t==2
        opti.set_initial(x, repmat(x_arr(1:n,t),N+1, 1));
        opti.set_initial(z, repmat(x_arr(1:n,t),N+1, 1));
        opti.set_initial(s, zeros(N+1,1));
        opti.set_initial(xr, x_arr(1:n,t));
        opti.set_initial(zr, x_arr(1:n,t));
        opti.set_initial(v, repmat(u_arr(1:m,t), N, 1));
        opti.set_initial(w, w_bar*ones(N,1));
    else
        opti.set_initial(x, [sol.value(x(n+1:end)); sol.value(x(end-2:end))]);
        opti.set_initial(z, [sol.value(z(n+1:end)); sol.value(z(end-2:end))]);
        opti.set_initial(u, [sol.value(u(m+1:end)); sol.value(u(end))]);
        opti.set_initial(v, [sol.value(v(m+1:end)); sol.value(v(end))]);
        opti.set_initial(xr, [sol.value(xr)]);
        opti.set_initial(zr, [sol.value(zr)]);
        opti.set_initial(w, [sol.value(w(2:end)); sol.value(w(end))]);
        opti.set_initial(s, [sol.value(s(2:end)); sol.value(s(end))]);
        opti.set_value(lambda_old,  opt_steady_online(n,m, theta_est, sol.value(u(end)), sol.value(xr), true, obj)  );
    end
    
    try
        tic;
        sol = opti.solve();
        Tsolve = toc
    catch ME
        % debugging information on constraints displayed
        disp("Error while solving iteration " + num2str(t));
        disp(ME)

        constraints = Function('constraints', {z,x,v,u,s,w, p, theta_LMS, xr, zr, ur, x_k}, {opti.debug.g});
        conval = full(constraints(opti.debug.value(z), opti.debug.value(x), opti.debug.value(v), opti.debug.value(u),...
            opti.debug.value(s), opti.debug.value(w), opti.debug.value(p), opti.debug.value(theta_LMS),...
            opti.debug.value(xr), opti.debug.value(zr), opti.debug.value(ur), opti.debug.value(x_k)))';
        [conval; concheck']
        objective = Function('objective', {x,xr,p}, {opti.debug.f});
        full(objective(opti.debug.value(x), opti.debug.value(xr), opti.debug.value(p)))
    end

    u_sol = sol.value(u);
    v_sol = sol.value(v);

    x_sol = sol.value(x);
    assert(norm(x_sol(1:n)-x_arr(:,t))<=1e-3)       % check validity of solution

    % true system propagation
    noise_true = uniform(W.A, W.b);                 % get noise
    x_next = syst.dynamic(x_arr(:,t), u_sol(1:m,1), h, theta_true) + noise_true;
        % update variables / storage
    z_arr(1:n, t+1) = syst.dynamic(z_arr(1:n,t), v_sol(1:m,1), h, theta_nom);
    x_arr(1:n, t+1) = x_next;
    u_arr(1:m, t) = u_sol(1:m);
    v_arr(1:m, t) = v_sol(1:m);
    %xr_arr(1:n, t+1) = sol.value(xr);
    %ur_arr(1:n, t+1) = sol.value(ur);
    %zr_arr(1:n, t+1) = sol.value(zr);
    %vr_arr(1:n, t+1) = sol.value(vr);
    w_arr(1:n, t)= noise_true;
    theta_arr(1:2, t+1) = theta_est;
    Tsolve_arr(:, t) = Tsolve;
    

    %%%% TIME VARYING PARAMETER UPDATE
    if TV 
       theta_true = theta_true+[0.00000025;-0.000004];
       LMS_params.const.theta_true = theta_true;
    end
end


%% SAVE results
switch type
    case "AEMPC"
        save(['Results/Simulations/AE_scaled--N_' , num2str(N),  '-Tsim_' , num2str(Tsim),  '-mu_', num2str(mu),  '-h_', num2str(h), '-LMS_', num2str(LMS), '-TV_', num2str(TV), '.mat'], "w_arr", "u_arr", "v_arr", "x_arr", "z_arr", "xr_arr", "zr_arr", "sol", "theta_arr", "h", "Tsim", "mu", "w_bar", "Lw", "P", "P_f", "x_min", "x_max", "Tsolve_arr");
    case "EMPC_estimate"
        save(['Results/Simulations/E_wrong--N_' , num2str(N),  '-Tsim_' , num2str(Tsim),  '-mu_', num2str(mu),  '-h_', num2str(h), '-LMS_', num2str(LMS), '-TV_', num2str(TV), '.mat'], "w_arr", "u_arr", "v_arr", "x_arr", "z_arr", "xr_arr", "zr_arr", "sol", "theta_arr", "h", "Tsim", "mu", "w_bar", "Lw", "P", "P_f", "x_min", "x_max", "Tsolve_arr");
    case "EMPC_perfect"
        save(['Results/Simulations/E_true--N_' , num2str(N),  '-Tsim_' , num2str(Tsim),  '-mu_', num2str(mu),  '-h_', num2str(h), '-LMS_', num2str(LMS), '-TV_', num2str(TV), '.mat'], "w_arr", "u_arr", "v_arr", "x_arr", "z_arr", "xr_arr", "zr_arr", "sol", "theta_arr", "h", "Tsim", "mu", "w_bar", "Lw", "P", "P_f", "x_min", "x_max", "Tsolve_arr");
end
disp(['INFO: Saved results of the simulation ', type])

end 


