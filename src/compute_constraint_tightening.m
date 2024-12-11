% compute_constraint_tightening.m


c = zeros(8,1);
c(1) = norm([1,0,0]*P^(-0.5),2);             c(4) = c(1);     % x_1
c(2) = norm([0,1,0]*P^(-0.5),2);             c(5) = c(2);     % x_2
c(3) = norm([0,0,1]*P^(-0.5),2);             c(6) = c(3);     % x_3
c(7:8) = norm(1*K*P^(-0.5),2)* ones(2,1);