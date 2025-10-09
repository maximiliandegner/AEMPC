% compute_terminal_cost.m
run parameter_def.m
theta_nom
load('opt_steady.mat', 'r_set')
load("Offline/grid_disc.mat", "P", "K", "h")
n = 3; m = 1;
z_s = full(r_set(1:3)); u_s = full(r_set(4));
x_b = [min(z_s)-0.1; max(z_s)+0.1];
[H,p_val] = syst.compute_maxH(z_s,u_s,theta_nom, h, x_b, K, {});

% \tilde{Q} = Q + M = Q + \sqrt{3}*H*\|p(\theta)\|
M = sqrt(n)*H*p_val;
assert(all(eig(M)>=0))

nonconst = false;
run constraint_def.m
run solve_terminalcost_SDP;

save("Offline/terminal_cost_P_f.mat", "P_f")                         
disp("End of terminal cost computation.")