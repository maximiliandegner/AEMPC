% LMI_improved.m

function [P, K, w_bar, rho, c_j,Lw] = LMI_improved(h, n_steps, RHO, LW_range, w_max, save_to_file, equilib, b_x, b_u)
%%LMI_IMPROVED compute the tube automatically for a range of RHO and Lw values
% The equilib input must be composed of the desired steady-state and
% steady-input, i.e., equilib = [x_s;u_s] \in \mathbb{R}^{n+m}.

n=3;m=1;  
mfilePath = mfilename('fullpath');
if ~isempty(mfilePath)
        mfilePath = mfilePath(1:find(mfilePath == '\', 2, 'last')-1);
end
addpath("\..\");
    % n_steps=10;%specify how many grid point for each dimension
    switch nargin
        case 0
            h = 0.025; %0.015
            n_steps = 3;
            RHO = linspace(0.97, 0.9999, 10); %0.985
            LW_range = linspace(0.001, 0.03, 10);
            w_max = 0.002; %005
            save_to_file = "";
            theta_true = syst.theta_true;
            theta_nom = syst.theta_nom;
            run([mfilePath, '\parameter_def.m'])
            run constraint_def.m
            run([mfilePath, '\opt_steady_offline.m']); load([mfilePath, '\opt_steady.mat'], "r_set"); equilib = full(r_set);
        case 5
            save_to_file = "";
            RHO = linspace(0.98, 0.9999, 6);
            theta_true = syst.theta_true;
            theta_nom = syst.theta_nom;
            run([mfilePath, '\parameter_def.m'])
            run constraint_def.m
            run([mfilePath, '\opt_steady_offline.m']); load([mfilePath, '\opt_steady.mat'], "r_set"); equilib = full(r_set);
        case 6
            disp("Saving to path\n"+save_to_file)
            theta_true = syst.theta_true;
            theta_nom = syst.theta_nom;
            run([mfilePath, '\parameter_def.m'])
            run constraint_def.m
            run([mfilePath, '\opt_steady_offline.m']); load([mfilePath, '\opt_steady.mat'], "r_set"); equilib = full(r_set);
        case 8    
            error("Wrong inputs: Please provide either 'b_x' and 'b_u' or neither.")
        case 9
            disp("Saving to path\n"+save_to_file)
        otherwise
            error("Wrong inputs.")
    end
    
    import casadi.*    
    a_1 = syst.scl(1); a_2 = syst.scl(2);      % scaling factors for parameter in model
    run([mfilePath, '\parameter_def.m'])
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

%% RHO loop
% initialize storage variables
P_arr = {};
K_arr = {};
input_tgt = {};
state_tgt = {};

for Lw = LW_range               % fixing a specific Lw value
    disp(['The next simulations use Lw = ', num2str(Lw)])
for rho_ii = 1:length(RHO)       % fixing a specific rho value
    rho = RHO(rho_ii);
    disp(['This solution uses rho = ', num2str(RHO(rho_ii))])

    if 1-rho-Lw < 0
        warning(['Lw value too large; rho=', num2str(rho), '; Lw=',num2str(Lw), '; 1-rho-Lw=', num2str(1-rho-Lw),'.']);
        pause(2);
        break;
    end

    % create the SDP and solve it 
    addpath("Offline\");
    run solve_offline_SDP;
    
    
    % obtain solution if solved successfully
    if solved
        Y=value(Y);  
        X=value(X);
        fprintf("Matrix P is:\n")
    
        if nonconst
            Y_0 = value(Y_0); Y_1=value(Y_1); Y_2=value(Y_2); Y_3=value(Y_3); Y_4=value(Y_4); Y_5=value(Y_5); Y_6=value(Y_6);  
            X_0 = value(X_0); X_1=value(X_1); X_2=value(X_2); X_3=value(X_3); X_4=value(X_4); X_5=value(X_5); X_6=value(X_6);
            X_max=value(X_max);
            P = value(inv(X_max))
        else
            P = value(inv(X))
        end
        K = value(Y*P)
        
        %store the output
        toc(t) 
        
        disp( ['Rho = ', num2str(RHO(rho_ii))] )
        disp("constraint-tightening terms")
        disp("input: " + num2str(norm(K*P^-0.5)))
        disp("state: " + num2str(norm(P^-0.5)))
        
        %save iterate solution
        P_arr{rho_ii} = P;                                  % P matrices
        P_maxeig{rho_ii} = max(eig(P));                     % max. eigenvalue of P matrices
        K_arr{rho_ii} = K;                                  % K matrices
        input_tgt{rho_ii} = norm(K*P^-0.5);                 % input tightenings
        state_tgt{rho_ii} = norm(P^-0.5);                   % state tightenings
        c_j_arr{rho_ii} = value(c_j);                       % values of variable c_j
        wbar_arr{rho_ii} = value(w_max);                    % values of w_max used
        Lw_arr{rho_ii} = value(Lw);                         % values of Lw used
        s_infty_arr{rho_ii} = value(w_max)/(1-rho-Lw);      % values of max. tube size
    
    else
        %save iterate solution if not solved successfully
        P_arr{rho_ii} = nan;
        P_maxeig{rho_ii} = nan;
        K_arr{rho_ii} = nan;
        input_tgt{rho_ii} = nan;
        state_tgt{rho_ii} = nan;
        c_j_arr{rho_ii} = nan;
        wbar_arr{rho_ii} = value(w_max);
        Lw_arr{rho_ii} = value(Lw);
        s_infty_arr{rho_ii} = value(w_max)/(1-rho-Lw);
    end %if solved
        

    end %for rho

input_str=input('Please hit ENTER to continue or enter *results* to store a solution:', 's');
if isempty(input_str)
    disp("Continuing with next Lw value...")

elseif contains(input_str, 'results')
    disp(['Maximum eigenvalues of the P matrices with LW = ', num2str(Lw)]); disp(P_maxeig);
    disp(1:length(cell2mat(P_maxeig)))
    it = input("Please enter the number of the solution that you want to save (0 to length(P_maxeig)):");
    P = P_arr{it}; K = K_arr{it}; rho = RHO(it);
    save("Offline/grid_disc.mat", 'Y_0','X_0','h', "rho", "K", "P", "w_max", 'Lw');
    w_bar = w_max;
    c_j = c_j_arr{it};
    return

else
    warning("Inadmissible input. Continuing with next Lw value...")
end

end %for Lw

it = 1;
P = P_arr{it}; K = K_arr{it}; rho = RHO(it);
save("Offline/grid_disc.mat", 'Y_0','X_0','h', "rho", "K", "P", "w_max", 'Lw');

% if save_to_file ~= ""
%     if contains(save_to_file, "cost")
%         save(save_to_file,'Y_0','Y_1','Y_2','Y_3','Y_4','Y_5','Y_6','X_0','X_1','X_2','X_3','X_4',...
%        'X_5','X_6','X_max','h','x_1max','x_1min','x_2max','x_2min','x_3max','x_3min',...
%        'u_max','u_min', "rho", "K", "P", "w_max")
%     else % e.g., tube offline design
%         w_max = sqrt(value(w_max));
%         disp("w_max is: " + num2str(value(w_max)))
%           save("Offline/" + save_to_file,'Y_0','Y_1','Y_2','Y_3','Y_4','Y_5','Y_6','X_0','X_1','X_2','X_3','X_4',...
%         'X_5','X_6','X_max','h','x_1max','x_1min','x_2max','x_2min','x_3max','x_3min',...
%         'u_max','u_min', "rho","K", "P", "w_max")
%     end
% end

% run('C:\Users\maxim\Documents\Masterthesis\adaptive-economic-mpc\matlab\nonlinear-pdep\fresh_trial\Simulations_CSTR_only_needed/opt_steady_Casadi.m')
end