% opt_steady_online.m
function [ell_out,obj] = opt_steady_online(n,m, theta, u_s, x_s, silent, obj)
import casadi.*
%constants
load("Offline\grid_disc.mat", "h")
alpha=2;delta=0.55;
a_1 = 10^5*theta(1); a_2 = 4*10^2*theta(2);

x_r1min=0.03;
x_r1max=0.8;
x_r2min=0.03;
x_r2max=0.8;
x_r3min=0.03;
x_r3max=0.8;
% fix the input
u_rmin=u_s;
u_rmax=u_s;
u_rmin = u_s;
u_rmax = u_s;

%% bounds
lb=[x_r1min;x_r2min;x_r3min;u_rmin];
ub=[x_r1max;x_r2max;x_r3max;u_rmax];

if ~isempty(obj.solver_lambda)
    solver = obj.solver_lambda;
else
    %%%%%% solver not yet initialized
    %define symbolic variables
    u = MX.sym('u',m,1);
    x = MX.sym('x',n,1);
    
    %stage cost
    cost = -x(2);  
    y=[x;u]; 
    %constraints= syst.fun(x,u,[theta(1);theta(2)]);
    constraints = h*fun(x,u, a_1, a_2, alpha, delta);
    %% NLP
    nlp = struct('x', y, 'f', cost, 'g', constraints);
    
    % Create IPOPT solver object
    if silent 
        opts_dict.ipopt.print_level = 0;
        opts_dict.print_time = 0;
    else
        opts_dict.ipopt.print_level = 5;
        opts_dict.print_time = 0;
    end
    % opts_dict.ipopt.tol = 1e-3;
    solver = nlpsol('solver', 'ipopt', nlp, opts_dict);
    obj.solver_lambda = solver;
    % Solve the NLP
end %if obj.solver_lambda ==

%%%%% solver opt problem
res = solver('x0' , [x_s;u_s],... % solution guess
             'lbx', lb,...           % lower bound on x
             'ubx', ub,...           % upper bound on x
             'lbg',    -1e-5*ones(3,1),...           % lower bound on g
             'ubg',    1e-5*ones(3,1));             % upper bound on g
%0);...%
 
%%
% r_set=res.x;
ell_out = full(-res.x(2));
obj.lambda = ell_out;

if ~contains(solver.stats().return_status, 'Solve_Succeeded')
    % if contains(solver.stats().return_status, 'Unset')
    %     load('opt_steady.mat', 'r_set')
    %     ell_out = full(-r_set(2));
    % end
    warning([solver.stats().return_status '     -- -- --     opt_steady_online.m']);
end


%% 
% save('opt_steady','r_set', 'l_set')
% toc(t)

function f=fun(x,u, a_1, a_2, alpha, delta)
% alpha=2;delta=0.55;% a_2=400;a_1=3*10^6;

f=[-a_1*exp(-1/x(3))*x(1)^alpha-a_2*exp(-delta/x(3))*x(1)-x(1)+1;...
    a_1*exp(-1/x(3))*x(1)^alpha-x(2);...
    -x(3)+u];
end 
end



