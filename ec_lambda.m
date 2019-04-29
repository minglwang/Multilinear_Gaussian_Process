function [p_std,p_ris]=ec_lambda(x,y)
fun_dc = @(p) p(1)*x./(p(2)+x).*1./(1+p(3)*x);
p0 = [1,1,1]; % intial values
%% standard quadratic approach
p_std = lsqnonlin(@(p)fun_dc(p)-y, p0, [-1e5,0,0]);
%% regularization approach with empirical choice of lambda
lambda=1e-3*mean(x);
 cost_G = @(p) norm(fun_dc(p)-y)^2+lambda*norm([p(2),p(3)],1);
 lb=[-1e5,0,0];
 p_ris = fmincon(cost_G,p0,[],[],[],[],lb);
end