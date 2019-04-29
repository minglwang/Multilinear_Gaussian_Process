function [best_para,best_lambda,lambda,valid_err]=cv_lambda(x,y,foldn,fun_dc)
% use a k-fold cross validation to find the penalty factor lambda as
% describe in Algorithm 4.
indices = crossvalind('Kfold',y,foldn);
lambda_range = (0.01:0.3:3);
p_g = zeros(length(lambda_range),3);
for j = 1:length(lambda_range)
lambda(j) = lambda_range(j);
test_err=0; para=zeros(1,3); cost_val=1e10; tol_fval=0;
for i = 1:foldn
    test = (indices == i); 
    train = ~test;
    x_train = x(train);
    x_test = x(test);
    ntest = length(x_test);
    % The penalized cost function in 
    cost_G = @(p) norm(fun_dc(p,x_train)-y(train))^2+lambda(j)*norm([p(2),p(3)],1);
    p0 = [1,1,1];
    [p_GP,fval] = fmincon(cost_G,p0,[],[],[],[],[-1e5,0,0]);
    if fval < cost_val
        cost_val = fval;
        para = p_GP;
    end
    test_err = test_err+norm((fun_dc(p_GP,x_test)-y(test)))^2;
end
    valid_err(j) = test_err; p_g(j,:) = para;
end
    [~,ind] = min(valid_err);
    best_para = p_g(ind,:); best_lambda = lambda(ind);
end