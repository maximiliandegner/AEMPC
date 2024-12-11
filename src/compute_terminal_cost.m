% compute_terminal_cost.m

%%%%%%%% Please add the folder "Offline" to the Matlab path and then run
%%%%%%%% this file from the nonlinear_sys folder.

run parameter_def.m
theta_nom
load('opt_steady.mat', 'r_set')
load("Offline/grid_disc.mat", "P", "K", "h")
n = 3;
z_s = full(r_set(1:3)); u_s = full(r_set(4));
x_b = [min(z_s)-0.1; max(z_s)+0.1];
% [A,B,X,Y,K] = syst.getXY(z_s,u_s,h, theta_nom); %P = inv(X);
% NOTE: Matlab's dlqr returns K such that u= - (Kx) is the best feedback and not +Kx!
P;
K;
[H,p_val] = syst.compute_maxH(z_s,u_s,theta_nom, h, x_b, K, {});

% \tilde{Q} = Q + M = Q + 3*H*\|p(\theta)\|
M = sqrt(n)*H*p_val;


assert(all(eig(M)>=0))
% Hessian of cost is zero.
[P_f, ~] = LMIs_grid_disc(h, zeros(3)+M , 0, 3, 0.999, "Offline/terminal_cost.mat");
                         % h, Q, R, n_steps, rho, save_to_file)

save("Offline/terminal_cost_P_f.mat", "P_f")                         
disp("End of terminal cost computation.")