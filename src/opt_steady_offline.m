%% setup
% t=tic;
import casadi.*
run parameter_def.m 

alpha = syst.alpha;
delta = syst.delta;
run parameter_def.m

load('Offline/grid_disc.mat', 'K')

x_r1min=0.05;
x_r1max=0.8;
x_r2min=0.05;
x_r2max=0.8;
x_r3min=0.05;
x_r3max=0.8;
u_rmin=0.049;
u_rmax=0.449;

%define symbolic variables
n = syst.n; m = syst.m;  % get system dimension
u = MX.sym('u',m);
x = MX.sym('x',n);
y=[x;u]; 

%stage cost
cost = -x(2);  
%dynamics
constraints=syst.fun(x,u,syst.theta_true);

%% bounds
lb=[x_r1min;x_r2min;x_r3min;u_rmin];
ub=[x_r1max;x_r2max;x_r3max;u_rmax];

%% NLP
nlp = struct('x', y, 'f', cost, 'g', constraints);

% Create IPOPT solver object
opts.ipopt.print_level = 0;
opts.print_time = 0;
solver = nlpsol('solver', 'ipopt', nlp, opts);

% Solve the NLP
res = solver('x0' , [0.35,0.35,0.35,0.35],...%[0.148498, 0.189011, 0.136954, 0.136958],... % solution guess
             'lbx', lb,...           % lower bound on x
             'ubx', ub,...           % upper bound on x
             'lbg',    -1e-5,...           % lower bound on g
             'ubg',    1e-5);             % upper bound on g
 
%% save results
r_set=res.x;
l_set = -res.x(2);

% toc(t)
v_set = res.x(4)-K*res.x(1:3);
save('opt_steady.mat','r_set', 'l_set', 'v_set')

