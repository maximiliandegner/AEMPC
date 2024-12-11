function [theta_hat, theta_nom] = param_estim(x, u, x_next, h, theta, uncert_set, params)
%PARAM_ESTIM Compute LMS and nominal parameter estimates

% Recall constants
mu = params.const.mu;
theta_true = params.const.theta_true;
p = length(params.const.theta_true);

% G_k = zeros(params.sys.n, p);
G_k = h*[-10^5*x(1)^2*exp(-1/x(3)), -4e2*x(1)*exp(-0.6/x(3));
          10^5*x(1)^2*exp(-1/x(3)),  0;
          0,                      0];
% G_k = G_k*diag([10,1]);


%% Point estimates 
%%% LMS update equations for best estimate

% LMS update equation
temp = x_next - (syst.dynamic(x,u,h,theta));
% temp2 = syst.dynamic(x,u,h,theta_true) - syst.dynamic(x,u,h,theta);
theta_tilde = theta + mu*G_k'*(temp);

% LMS projection
try
    ret = uncert_set.project(theta_tilde);
    theta_hat = ret.x; ok=ret.exitflag; dist=ret.dist;
    if ok == 0
        error('Projection of LMS estimate failed')
    end %if
catch ME
    disp(ME)
end

%%% Nominal parameter as center of hyperbox
%%% only for set-membership
% outBox = uncert_set.outerApprox;
% theta_nom = sum(outBox.A.*outBox.b)'/2;
theta_nom = params.const.theta_nom;

end %function