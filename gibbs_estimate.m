function [v_hat, theta] = gibbs_estimate(y,x) 
  
    % algorithmetic parameters
    EM_steps = 100; % Number of EM updates
    rel_tol = 1e-1; % Relative tolerance of the updates
    Burnin1 = 20; % MC parameters
    MC_steps1 = 100;
    Burnin2 = 50;
    MC_steps2 = 1000;
    N = length(y);
    m=size(x,2);
    H=zeros(N,N,m);
    Kernel=zeros(N,N,m);
    S = zeros(N,N,m);
    f=rand(N,m);
    vi=zeros(N,1);
    
    % Initial values
    ini_state = ones(N,1);
    err = 0;
    sigma2 = .01;
    lambda = ones(m,1)*0.1;%
    rho = ones(m,1)*0.05;
    theta = [lambda;rho;sigma2];
    
    % kernel design
    K = @(x,y,r) r(1)*(x(:)./y(:)').^(-r(2)*log(x(:)./y(:)')); % kernel in (26)

    for kk=1:m
    H(:,:,kk)= K(x(:,kk), x(:,kk), [1,1]);
    Kernel(:,:,kk)=lambda(kk)*H(:,:,kk).^rho(kk);
    end
    
    disp('Starting EM iterations');
    disp('---------------');
    disp('| iter: tol   |')
    disp('---------------');
    for jj = 1:EM_steps
        
        gamma_jj = 1/(jj^.6);
        theta_old = theta;
        
        % E-STEP
        for kk=1:m
        Kernel(:,:,kk)=lambda(kk)*H(:,:,kk).^rho(kk);
        end
        S_step = zeros(N,N,m);
        err_step = 0;
        
        for mm = -Burnin1:MC_steps1
            for kk=1:m
                condition_set=setdiff(1:m,kk);
                f_cond_on=prod(f(:,condition_set),2);
                vi=sample(Kernel(:,:,kk), y, f_cond_on, sigma2);
                f(:,kk)=vi;
            end
            if mm >0
                Burnin1 = 0;
                for kk=1:m
                S_step(:,:,kk) = S_step(:,:,kk) + f(:,kk)*f(:,kk)'/MC_steps1;
                end
                err_step = err_step + norm(y-prod(f,2))^2/MC_steps1;
            end
        end
        for kk=1:m
        S(:,:,kk) = (1-gamma_jj)*S(:,:,kk) + gamma_jj*S_step(:,:,kk);
        end
        err = (1-gamma_jj)*err + gamma_jj*err_step;
        
        % M-step: update hyperparameters
        for kk=1:m
        [lambda_step, rho_step] = optimizeHypers(S(:,:,kk), H(:,:,kk));
        lambda(kk)=lambda_step;
        rho(kk)=rho_step;
        end
        sigma2 = err/N;
        
        theta = [lambda; rho; sigma2];
        tol =norm(theta - theta_old) / norm(theta_old);

        disp(sprintf('| %4d: %0.3f |',jj,tol));
        if tol < rel_tol
            reason = 'Convergence: Tolerance met!';
            break
        end
    end
    if ~exist('reason')
        reason = 'Maximum iterations reached!';
    end
    disp(sprintf(['---------------\n', reason]));
    for kk=1:m
        Kernel(:,:,kk)=lambda(kk)*H(:,:,kk).^rho(kk);
    end
    v_hat = zeros(N,m);
    for mm = -Burnin2:MC_steps2
        for kk=1:m
                condition_set=setdiff(1:m,kk);
                f_cond_on=prod(f(:,condition_set),2);
                vi=sample(Kernel(:,:,kk), y, f_cond_on, sigma2);
                f(:,kk)=vi;
        end
        if mm >0
            v_hat = v_hat + f/MC_steps2;
        end
        
    end
end

function [lambda_out, rho_out] = optimizeHypers(S, H)

    rho_grid = linspace(0.0001,1.,300);
    best_cost = Inf;
    N = size(S,1);


    for rho = rho_grid
        [U, D] = svd(H.^rho);
        D = diag(D);
        n = find(sum(D) - cumsum(D) < 1e-9,1);
        tr = trace( diag(1./D(1:n))*U(:,1:n)'*S*U(:,1:n));
        ld = sum(log(D(1:n)));

        cost = ld + N*log(tr);
        if cost < best_cost
            best_cost = cost;
            rho_out = rho;
            lambda_out = tr/N;
        end

    end
end

function v = sample(K, y, m, sigma2)
    N = length(y);
    P = (K*diag(m.^2)/sigma2 + eye(N))\K;
    P = .5*(P + P');
    m = P * (m.*y)/sigma2;
    v = m + chol(P+1e-6*eye(N))'*randn(N,1);
end

