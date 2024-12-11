% syst.m

classdef syst
    methods (Static)
        function a = static_param()
            d = 2;
            a = zeros(d,1);
            a(1) = 1;% 10^5; --> scaling of 10^5 included in model equations
            a(2) = 1; % 400; --> scaling of 4*10^2 included in model equations
        end

        function [A,B]=getA_d(x,u,h, theta)
            import casadi.*
            % k1=fun(x,u);
            %RuKu1 (Euler) discretization
            %1. derivative k_1: [A_c,B_c]
            switch nargin
                case 3
                    [A_1,B_1]=syst.getA_c(x,u);
                case 4
                    [A_1,B_1]=syst.getA_c(x,u, theta);
                otherwise
                   error("nargin not right.")
            end
            n=3;

            % assemble A and B matrices
            A=eye(n)+h*A_1;
            B=h*B_1;
        end %function
            
        function [A,B]=getA_d_RK4(x,u,h, theta)
            import casadi.*
            %Jacobian of discrete-time system (hard coded using cont.-time Jacobian +
            %RK4)
            n=3; 
            switch nargin
                case 3
                    k1=syst.fun(x,u);
                    k2=syst.fun(x+h/2*k1,u);
                    k3=syst.fun(x+h/2*k2,u);
                    k4=syst.fun(x+k3,u);
                    %RuKu4 discretization
                    %1. derivative k_1: [A_c,B_c]
                    [A_1,B_1]=syst.getA_c(x,u);
                    %2. derivative k_2:
                    [A_2,B_2]=syst.getA_c(x+h/2*k1,u);
                    A_k2=A_2*(eye(n)+h/2*A_1);
                    B_k2=B_2+h/2*A_2*B_1;
                    %3. derivative k_3:
                    [A_3,B_3]=syst.getA_c(x+h/2*k2,u);
                    A_k3=A_3*(eye(n)+h/2*A_k2);
                    B_k3=B_3+h/2*A_3*B_k2;
                    %4. derivative k_4:
                    [A_4,B_4]=syst.getA_c(x+h*k3,u);
                    A_k4=A_4*(eye(n)+h*A_k3);
                    B_k4=B_4+h*A_4*B_k3;
                case 4
                    k1=syst.fun(x,u,theta);
                    k2=syst.fun(x+h/2*k1,u,theta);
                    k3=syst.fun(x+h/2*k2,u,theta);
                    k4=syst.fun(x+k3,u,theta);
                    %RuKu4 discretization
                    %1. derivative k_1: [A_c,B_c]
                    [A_1,B_1]=syst.getA_c(x,u,theta);
                    %2. derivative k_2:
                    [A_2,B_2]=syst.getA_c(x+h/2*k1,u,theta);
                    A_k2=A_2*(eye(n)+h/2*A_1);
                    B_k2=B_2+h/2*A_2*B_1;
                    %3. derivative k_3:
                    [A_3,B_3]=syst.getA_c(x+h/2*k2,u,theta);
                    A_k3=A_3*(eye(n)+h/2*A_k2);
                    B_k3=B_3+h/2*A_3*B_k2;
                    %4. derivative k_4:
                    [A_4,B_4]=syst.getA_c(x+h*k3,u,theta);
                    A_k4=A_4*(eye(n)+h*A_k3);
                    B_k4=B_4+h*A_4*B_k3;
                otherwise
                    error("nargin mismatch.")
            end
            
            A=eye(n)+h/6*(A_1+2*A_k2+2*A_k3+A_k4);
            B=h/6*(B_1+2*B_k2+2*B_k3+B_k4);
        end
            
        function [A,B]=getA_c(x,u, theta)
            import casadi.*
            switch nargin
                case 1
                    error("Not enough arguments.")
                case 2
                    a = syst.static_param();
                    a_1 = a(1); a_2 = a(2);
                case 3
                    a_2=theta(2); a_1=theta(1);
            end
            alpha=2;  delta=0.55;
            %including scaling
            a_1 = a_1*10^5;
            a_2 = a_2*4*10^2;

            %Jacobian of cont.-time system
            n=3;
            m=1;
            %compute lineariztion
            A=zeros(n);B=zeros(n,m);
            if ~isa(x,'double')
                A = MX(A); B = MX(B);
            end

            A(1,1)=-1-alpha*a_1*x(1)^(alpha-1)*exp(-1/x(3))-a_2*exp(-delta/x(3));
            A(1,2)=0;
            A(1,3)=-a_2*x(1)*delta/(x(3)^2)*exp(-delta/x(3))-a_1*x(1)^alpha/(x(3)^2)*exp(-1/x(3));

            A(2,1)=alpha*a_1*x(1)^(alpha-1)*exp(-1/x(3));
            A(2,2)=-1;
            A(2,3)=a_1*x(1)^alpha/(x(3)^2)*exp(-1/x(3));

            A(3,1)=0;
            A(3,2)=0;
            A(3,3)=-1;
            
            B=[0;0;1];
        end
            
        function f=fun(x,u,theta)
            import casadi.*
            switch nargin
                case 1
                    error("Not enough arguments.")
                case 2
                    a = syst.static_param();
                    a_1 = a(1); a_2 = a(2);
                case 3
                    a_2=theta(2); a_1=theta(1);
            end
            alpha=2;  delta=0.55;
            % including scaling
            a_1 = 10^5*a_1;
            a_2 = 4*10^2*a_2;
            %continuous-time dynamics
            f=zeros(size(x,1),1);
            if isa(x,"MX")
                f=MX(f);
            elseif isa(x,"SX")
                f=SX(f);
            end
            f(1)=-a_1*exp(-1/x(3))*x(1)^alpha-a_2*exp(-delta/x(3))*x(1)-x(1)+1;
            f(2)=a_1*exp(-1/x(3))*x(1)^alpha-x(2);
            f(3)=-x(3)+u;
        end 
            
        function x_new=dynamic(x,u,h, theta)
            import casadi.*
            %nonlinear discrete-time model
            %here, obtain through RK4 of cont.-time
            switch nargin
                case 3
                    k1=syst.fun(x,u);
                case 4
                    k1=syst.fun(x,u,theta);
                otherwise
                    error("nargin mismatch.")
            end
            x_new=x+h*k1;
        end
            
        function x_new=dynamic_RK4(x,u,h, theta)
            import casadi.*
            %nonlinear discrete-time model
            %here, obtain through RK4 of cont.-time
            switch nargin
                case 3
                    k1=syst.fun(x,u);
                    k2=syst.fun(x+h/2*k1,u);
                    k3=syst.fun(x+h/2*k2,u);
                    k4=syst.fun(x+h*k3,u);
                case 4
                    k1=syst.fun(x,u,theta);
                    k2=syst.fun(x+h/2*k1,u,theta);
                    k3=syst.fun(x+h/2*k2,u,theta);
                    k4=syst.fun(x+h*k3,u);
                otherwise
                    error("nargin mismatch.")
            end

            x_new=x+h/6*(k1+2*k2+2*k3+k4);
        end

        function [A,B,X,Y,K] = getXY(x_r,u_r,h, theta)
            import casadi.*
            switch nargin
                case 3
                    load Offline/grid_disc
                    [A,B]  = syst.getA_d(x_r,u_r,h);
                case 4
                    load Offline/grid_disc
                    [A,B]  = syst.getA_d(x_r,u_r,h, theta);
                otherwise
                    error("nargin mismatch.")
            end
            
            t1=A(1,1);t2=A(1,3);t3=A(2,1);t4=A(2,3);t5=B(1);t6=B(2);
             %calculate Y and X for every point depending on the parameters at that point
             Y = Y_0;%+Y_1*t1+Y_2*t2+Y_3*t3+Y_4*t4;%+Y_5*t5+Y_6*t6;
             X = X_0;%+X_1*t1+X_2*t2+X_3*t3+X_4*t4;%+X_5*t5+X_6*t6;
             
             K=Y*inv(X);
        end

        function Lw = compute_Lw(Theta, theta_bar, W, P, X_poly, h)
            %compute_Lw computes the Lipschitz constant L_w by gridding / brute force
            % Lipschitz constant for the uncertainty bound computed offline, given 
            % the terminal set defining matrix P, the parametric uncertainty set 
            % Theta and the additive disturbance set W
            import casadi.*
            P_sqrt = (P^(0.5)); P_isqrt = (inv(P_sqrt));
            d = size(Theta.A,2); n = size(W.A,2);
            num = 10;
            
            x_space1 = linspace(-X_poly.b(n+1), X_poly.b(1), num);
            x_space2 = linspace(-X_poly.b(n+2), X_poly.b(2), num);
            x_space3 = linspace(-X_poly.b(n+3), X_poly.b(3), num);

            %%%% For manual sanity checks
            % Lw_x = zeros(3,1); Lw_t = zeros(5,1);
            
            % setup the Lw computation
            Lw = zeros(num,1);
            j1max = size(Theta.V,1);
            j2max = size(W.V,1);
            Theta_vert = Theta.V;
            W_vert = W.V;
            x = SX.sym('x',3,1); t = SX.sym('t',5,1); 
            G= SX.zeros(3,2);
            G(1,1)=h*(-1e4*exp(-1/x(3))*x(1)^2);
            G(2,1)=h* (   1e4*exp(-1/x(3))*x(1)^2  );
            G(3,1)=0;
            G(1,2) = h*  ( -4*1e2*exp(-0.55/x(3))*x(1)  );
            G(2,2) = 0;
            G(3,2) = 0;
            fp= Function('fp', {x,t}, {P_sqrt*(  jacobian(G(:,1),x)*P_isqrt*t(1)+ jacobian(G(:,2),x)*P_isqrt*t(2)  )},...
                        {'x','t'}, {'o'});
            fp_test = Function('fp_test', {x}, {jacobian(G(:,1),x)}, {'x'}, {'o'});
            % loop over the whole x_space and select the highest Lipschitz constant
            for i = 1:num
                for j = 1:num
                    for k = 1:num
                        for j1=1:j1max
                            for j2 = 1:j2max

                                t_curr = [Theta_vert(j1,:)'-theta_bar; W_vert(j2,:)'];
                                temp = norm(full(fp([x_space1(i);x_space2(j);x_space3(k)], t_curr)));
                                
                                if temp > Lw(i)
                                    Lw(i) = temp;
                                    %%%% for manual sanity checks
                                    if temp >= 0.5
                                        Lw_x = [x_space1(i);x_space2(j);x_space3(k)]
                                        Lw_t = t_curr
                                        fp_test(Lw_x)
                                    end
                                end %if

                            end%j2
                        end%j1
                    end%for k
                end%for j
            end%for i
            Lw = max(Lw);
        % return Lw
    end %function compute_Lw
        
        function [H_hat,p_val] = compute_maxH(x_val,u_val,theta, h, x_bounds, F, obj)
            import casadi.*
            n = size(x_val,1);
            m = size(u_val,1);
            x = SX.sym('x',n,1);
            u = SX.sym('u',m,1);
            par = SX.sym('par', size(theta,1),1);
            
            dyn = syst.dynamic(x,u,h,par);

            H1 = Function('H1', {x,u,par}, {0.5*hessian(dyn(1),x)}, {'x','u', 'par'}, {'o'}); % *(x-x_val)*1/(2*1)
            H2 = Function('H2', {x,u,par}, {0.5*hessian(dyn(2),x)}, {'x','u', 'par'}, {'o'});
            H3 = Function('H3', {x,u,par}, {0.5*hessian(dyn(3),x)}, {'x','u', 'par'}, {'o'});

            num = 5;
            xr_1min=0.03;
            xr_1max=0.2;    % 0.4
            xr_2min=0.03;
            xr_2max=0.2;    % 0.2
            xr_3min=0.03;
            xr_3max=0.2;    % 0.2
            ur_min=0.049;
            ur_max=0.449;

            x_space1 = linspace(xr_1min, xr_1max, num);
            x_space2 = linspace(xr_2min, xr_2max, num);
            x_space3 = linspace(xr_3min, xr_3max, num);
            p_val = 0; m_val = 0;
            
            run parameter_def.m

            for i = 1:num
            % parfor i = 1:num
                for j = 1:num
                    for k = 1:num
                        for iii = 1:3
                            % for jjj = 1:3
                                theta_arr = (Theta_0.V(iii)'-theta_nom).*[1e4;4e2];
                                xi = [x_space1(i);x_space2(j);x_space3(k)];%-x_val;
                                ui = 0;
                                M1 = full(H1(xi,ui,Theta_0.V(iii)'));
                                M2 = full(H2(xi,ui,Theta_0.V(iii)'));
                                M3 = full(H3(xi,ui,Theta_0.V(iii)'));
                                m1 = max(eig(M1));
                                m2 = max(eig(M2));
                                m3 = max(eig(M3));
                                if max([m1,m2,m3])/2 >= m_val
                                    m_val = max([m1,m2,m3])/2;
                                end

                        
                                [A,B] = syst.getXY(xi,ui,h, theta_arr);
                                % disp(max(eig(A+B*F)))
                                temp_M = inv(eye(3)-(A+B*F));
                                pt = norm([0,-1,0]*temp_M);
                                if pt >= p_val
                                    p_val = pt;
                                end %if
                            % end %for jjj
                        end %for iii
                    end %for k
                end %for j
            end %for i
        disp(m_val)
        H_hat = m_val*eye(3);
        end

        function t = max_vert(Theta,W,nom)
            % vertices of Theta and W
            t1 = [Theta.V-nom';
                    zeros(size(W.V,1)-size(Theta.V,1), size(Theta.V,2))];
            t2 = W.V;
            t = [t1, t2]; t = max(t,[], 1);
        end %function max_vert

        function d = max_dist(Poly, M)
            %max_dist computes the maximum distance within a given polytop
            %Poly with the square of the matrix-weighted norm.
            d = 0;
            for i = 1:size(Poly.V,1)
                for j = i:size(Poly.V,1)
                    a = Poly.V(i, :); b = Poly.V(j, :);  
                    curr = (a-b)*M*(a-b)';
                    if d <= curr
                        d = curr;
                    end
                end
            end
        end% function max_dist

    end %method(static)
end %class

