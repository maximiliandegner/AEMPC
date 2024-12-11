%%% this file creates the offline SDP for computing P and K, and solves it    

%% definition of the decision variables
Y_0 = sdpvar(m,n);
X_0 = sdpvar(n);
if nonconst % if nonconstant matrices are desired
    Y_1=sdpvar(m,n); Y_2=sdpvar(m,n); Y_3=sdpvar(m,n); Y_4=sdpvar(m,n); Y_5=sdpvar(m,n);Y_6=sdpvar(m,n);
    X_1 = sdpvar(n);X_2 = sdpvar(n);X_3 = sdpvar(n);X_4 = sdpvar(n);X_5 = sdpvar(n);X_6 = sdpvar(n);
end
X_max=sdpvar(n);
c_j = sdpvar(length(h_xu),1);

% constraints 
con=[];
t=tic;

%% looping over all states and inputs
ii=0;%this counter shows progress to give you an idea how long it will take 
x_1max = b_x(1); x_1min = b_x(n+1);
x_2max = b_x(2); x_2min = b_x(n+2);
x_3max = b_x(3); x_3min = b_x(n+3);
u_max = b_u(1); u_min = b_u(m+1);

for vertT = 1:size(Theta_0.V,1)
theta = Theta_0.V(vertT,:);
theta_ = (theta - theta_nom') .*[a_1, a_2];      % compute distance from nominal parameter (e.g., center) to vertices and scale according to model

 for x1 = linspace(x_1min,x_1max,n_steps)
 for x2 = (x_2max+x_2min)/2 %hard coded here that x_2 doesn't matter
 for x3 = linspace(x_3min,x_3max,n_steps)
 for u =  linspace(u_min,u_max,n_steps)
     x=[x1;x2;x3];
     xplus = syst.dynamic(x,u,h,theta);%grid state-input-> simulate one step
     if xplus(1)>=x_1min&&xplus(3)>=x_3min&& xplus(1)<=x_1max&&xplus(3)<=x_3max
     
        %only look at (x,u), where also next state lies in constraints
        [A,B]  = syst.getA_d(x,u,h, theta);%compute Jacobian
    
         %calculate Y and X for every point depending on the parameters at that point
         Y = Y_0;
         X = X_0;
         if nonconst
            t1=A(1,1);t2=A(1,3);t3=A(2,1);t4=A(2,3);t5=B(1);t6=B(2);%syst.get nonlinear entries from (A,B) to parametrize P,K
            Y = Y+Y_1*t1+Y_2*t2+Y_3*t3+Y_4*t4+Y_5*t5+Y_6*t6;
            X = X++X_1*t1+X_2*t2+X_3*t3+X_4*t4+X_5*t5+X_6*t6;
            con=[con;X_max<=X];
         end

         for uplus =  linspace(u_min,u_max,2)
             xplusplus=syst.dynamic(xplus,uplus,h, theta);
             if xplusplus(1)>=x_1min&&xplusplus(3)>=x_3min&& xplusplus(1)<=x_1max&&xplusplus(3)<=x_3max
                 ii=ii+1;
    
                 [Ap,Bp]  = syst.getA_d(xplus,uplus,h, theta);
                 Xp = X_0;
             
                 if nonconst
                    tp1=Ap(1,1);tp2=Ap(1,3);tp3=Ap(2,1);tp4=Ap(2,3);tp5=Bp(1);tp6=Bp(2);
                    Xp = Xp+X_1*tp1+X_2*tp2+X_3*tp3+X_4*tp4+X_5*tp5+X_6*tp6;
                 end

                 ineq1=[rho^2*X,     (A*X+B*Y)';... %,     Q_eps*X;...
                        A*X+B*Y,        Xp];%,        zeros(n,n);... 
                        %X*Q_eps,    zeros(n),       eye(n)];
         
                 G_th = full(G_1x(x))*theta_(1) + full(G_2x(x))*theta_(2);
                 ineq11 = [Lw^2*X,        G_th*X;
                          (G_th*X)',     X];
        
                 con=[con;ineq1>=0];
                 con=[con;ineq11>=0];
             end %if xplusplus(1)

         end %for uplus
     end %if xplus
 end %for u
     
 end %for x3
 end %for x2
 end %for x1

% require the uncertainty bound only at the desired steady-state from 'equilib' 
for vertW = 1:size(W.V, 1)
     G_fval = full(G_f(equilib(1:n)));
     ineq2 = [X,    G_fval*theta_'+h*W.V(vertW,:)';
              (G_fval*theta_'+h*W.V(vertW,:)')' , w_max^2];
     con = [con;ineq2>=0];
end %for vertW

end %for vertT

% require the constraint tightening only at the desired equilibrium 'equilib'
for jj = 1:numel(h_xu)
con = [con; c_j(jj)/(-g_val(jj))*( w_max/(1-rho-Lw) ) <= 1];

ineq3 = [c_j(jj)           ,  H_xu(jj,:)*[X;Y];
         (H_xu(jj,:)*[X;Y])' ,  X];
con = [con; ineq3 >= 0; c_j(jj) >= 0];
end %for jj
disp("Number of constraints added: " + num2str(ii+jj))

% adding PSD of P 
con=[con;   X<= eye(n)];

disp('Finished adding constraints.')

%% running optimization, looping over RHO
disp(['Starting optimization for rho= ' num2str(rho)]);
options = sdpsettings('solver','mosek', 'verbose', 0);
cost = 0;

 weight = ones(8,1);%10*ones(length(c_j),1); weight(7:8) = [1;1];
 for jj=1:numel(c_j)
     cost = cost + weight(jj)*c_j(jj)^2*w_max^2/(1-rho-Lw)^2 * 1/(-g_val(jj))^2; 
 end

% % for debugging unbounded issues with YALMIP
con = [con; -10^8 <= allvariables(cost) <=10^8];    

ret = optimize(con,cost, options);

% throw warning if problem not solved
switch ret.problem
    case 0
        disp("Solved to optimality.")
        solved = true;
    otherwise
        warning(ret.info)
        solved = false;
end

