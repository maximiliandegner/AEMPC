% uncertainty_fcn.m

function w = uncertainty_fcn(x, u, h, theta, Theta, W, P)
    Psq = P^0.5;
    w = 0;
    delta = 0.55;
    a_1 = 10^5; a_2 = 4*10^2;
    p  = size(Theta.V, 2);
    nu = size(W.V,2);
    
    import casadi.*
    xx = MX.sym('xx', 3,1);
    G = h* [-xx(1)^2*exp(-1/xx(3)), -xx(1)*exp(-delta/xx(3));
        xx(1)^2*exp(-1/xx(3))  ,  MX(0);
        MX(0)                  ,  MX(0)];
    G_f = Function('G_f', {xx}, {G});
    A = [Theta.A, zeros(size(Theta.A,1),nu); zeros(size(W.A,1), p), W.A];
    b = [Theta.b; W.b];
    
    temp = [];
    for t1 = 1:size(Theta.V,1)
        param = (Theta.V(t1,:)'-theta) .* [a_1; a_2];
        for t2 = 1:size(W.V,1)
            noise = h*W.V(t2,:)';
            temp = [temp; sqrt( (G_f(x)*param + noise)'*P*(G_f(x)*param + noise) )];
            % if temp > w
            %     w = temp;
            % end
        end
    end
    w = temp;
end

