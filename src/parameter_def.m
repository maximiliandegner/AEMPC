% state and input constraints
x_min=0.03;
x_max=0.25;
u_min=0.049;
u_max=0.449;


% theta_nom  = syst.static_param.*[1.00; 1];% 
% theta_nom  = syst.static_param;
% d = length(theta_true);
% theta_true = [1.0049;0.9959].*theta_true;

% parameter polytope for projection in the LMS
A_theta = [eye(2); -eye(2)];
b_theta = [theta_nom(1)*1.015; theta_nom(2)*1.015; -theta_nom(1)*0.985; -theta_nom(2)*0.985];
% b_theta = [theta_true(1); theta_true(2); -theta_true(1); -theta_true(2)];
Theta_0 = Polyhedron(A_theta, b_theta);

assert(Theta_0.contains(theta_true))

% disturbance polytope
A_w = [eye(3); -eye(3)];
b = [0.0005; 0.0005; 0.0005];
b_w = [b; b];% [0.005*one(3,1); 0.005*ones(3,1)];
W = Polyhedron(A_w,b_w);

