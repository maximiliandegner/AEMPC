% t=tic;
import casadi.*
run parameter_def.m
%constants
alpha=2;delta=0.55;
global a_2 
global a_1
a_2=400;a_1=10^5;
run parameter_def.m
a_1 = theta_true(1)*a_1;
a_2 = theta_true(2)*a_2;


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
n=3;m=1;
u = MX.sym('u',m);
x = MX.sym('x',n);
%stage cost
cost = -x(2);  

y=[x;u]; 
constraints=fun(x,u);
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
 
%%
r_set=res.x
l_set = -res.x(2);

%% 
% toc(t)
v_set = res.x(4)-K*res.x(1:3);
disp(v_set)
save('opt_steady.mat','r_set', 'l_set', 'v_set')


function f=fun(x,u)
alpha=2;delta=0.55;% a_2=400;a_1=3*10^3;
global a_2 a_1
f=zeros(size(x,1),1);
f=[-a_1*exp(-1/x(3))*x(1)^alpha-a_2*exp(-delta/x(3))*x(1)-x(1)+1;...
    a_1*exp(-1/x(3))*x(1)^alpha-x(2);...
    -x(3)]...
    + [0;0;1]*u;
% f=[-a_1*exp(-1*x(3))*x(1)^alpha-a_2*exp(-delta*x(3))*x(1)-x(1)+1;...
%     a_1*exp(-1*x(3))*x(1)^alpha-x(2);...
%     -1/x(3)+u];
end 
