theta_true = syst.static_param();
theta_true(1) = 0.1*theta_true(1);
theta_nom  = theta_true.*[0.995; 1.005];% [1.005; 0.995];

% theta_nom  = syst.static_param.*[1.00; 1];% 
% theta_nom  = syst.static_param;
d = length(theta_true);
% theta_true = [1.0049;0.9959].*theta_true;

A_theta = [eye(2); -eye(2)];
b_theta = [theta_nom(1)*1.015; theta_nom(2)*1.015; -theta_nom(1)*0.985; -theta_nom(2)*0.985];
% b_theta = [theta_true(1); theta_true(2); -theta_true(1); -theta_true(2)];
Theta_0 = Polyhedron(A_theta, b_theta);

assert(Theta_0.contains(theta_true))

A_w = [eye(3); -eye(3)];
b = [0.0005; 0.0005; 0.0005];
b_w = [b; b];% [0.005*one(3,1); 0.005*ones(3,1)];
W = Polyhedron(A_w,b_w);

