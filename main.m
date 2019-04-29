%% main code for Algorithm 3
% {v_max, k_act, k_inh} are the model parameters;
% {h_act, h_inh, fun_dc} are the modulation functions;
% {F_v} is the true model;
% {N} smaple number;
% {epsilon} is a small positive value close to 0;
% {f_hat, g_hat} are the GP estimates.
% {y_fit} is the goodness of fit;
% {pf_std, pg_std} are the parameters identified using standard quadratic cost; 
% {pg_ris_cv, pg_ris_cv, pg_ris_ec, pf_ris_ec} are the parameters
% identified by L1-RIS via cross validation and empirical choice;
% {F, G, V} are the selected functions
% {p0} is the initial values
% {p_ft} is the fine-tuned parameters 

close all; clear all; clc
v_max = 1;
k_act = 12;
k_inh = 1/56; % true values k_inh=1/k2;
h_inh = @(x,k) 1./(1+k*x); 
h_act = @(x,k) x./(k+x); 
fun_dc = @(p,x) p(1)*x./(p(2)+x).*1./(1+p(3)*x);
F_v = @(c1,c2,k) v_max*h_act(c1,k(1)).*h_inh(c2,k(2));

N = 100;
c1min = 1;
c2min = 1;
c1max = 100;
c2max = 300;
c1 = c1min+c1max*rand(N,1);
c2 = c1min+c2max*rand(N,1);
epsilon=1e-3; 
v = F_v(c1,c2,[k_act,k_inh]);
sigma2 = 0.1*var(v);
y = v + sqrt(sigma2)*randn(N,1);

%% Step 1: multilinear GP
x=[c1,c2];
[v_hat, theta] = gibbs_estimate(y,x);
f_hat = v_hat(:,1); g_hat=v_hat(:,2);
y_hat = f_hat.*g_hat;
y_fit = 1-norm(y_hat - v)/norm(v - mean(v));

%% Step 2: Selection
% via 3-fold cross validation
[pf_ris_cv, lamda_f]=cv_lambda(c1,f_hat,3,fun_dc); [pg_ris_cv, lambda_g]=cv_lambda(c2,g_hat,3,fun_dc); 

% via empirical choice
[pf_std,pf_ris_ec]=ec_lambda(c1,f_hat); [pg_std,pg_ris_ec]=ec_lambda(c2,g_hat);

% model selection
H_act(1) = double(pf_ris_ec(2)>epsilon); H_inh(1) = double(pf_ris_ec(3)>epsilon);
H_act(2) = double(pg_ris_ec(2)>epsilon); H_inh(2) = double(pg_ris_ec(3)>epsilon);

ind_f=(H_act(1)==1) && (H_inh(1)==0); % 1 for true, 0 for false
ind_g=(H_act(2)==0) && (H_inh(2)==1); 

%% Step 3: Fine-tuning
F = @(k) (h_act(c1,k(1))).^(H_act(1)).*(h_inh(c1,k(2)).^(H_inh(1)));
G = @(k) ((h_act(c2,k(1))).^(H_act(2))).*(h_inh(c2,k(2)).^(H_inh(2)));
V = @(p) p(1)*F([p(2),p(3)]).*G([p(4),p(5)]);
p0 = [pf_ris_ec(1)*pg_ris_ec(1),pf_ris_ec(2:end),pg_ris_ec(2:end)];
lb = zeros(size(p0));
p_ft = lsqnonlin(@(p) V(p)-y,p0,lb);


%% Direct identification (for comparison)
V_d = @(p) p(5)*h_act(c1,p(1)).*h_inh(c1,p(2)).*h_act(c2,p(3)).*h_inh(c2,p(4));
p_d = lsqnonlin(@(p) V_d(p)-y,[1,1,1,1,1],[0,0,0,0,0]);

%% show the GP results
c1test=linspace(c1min+3, c1max, 20);
c2test=linspace(c2min+3, c2max, 20);
f_hat_test = interp1(c1, f_hat, c1test, 'pchip','extrap');
g_hat_test = interp1(c2, g_hat, c2test, 'pchip','extrap');
figure(1)
subplot(2,2,1)
plot(c1test,h_act(c1test,k_act),'r-','linewidth',1.5)
hold on
plot(c1test, f_hat_test/pf_ris_ec(1),'bo','linewidth',1)
% ylim([0,1])
xlabel('$c_1$','interpreter','latex')
ylabel('$h_1(c_1)$','interpreter','latex')
legend({'true','GP'},'location','southeast')
title('(a)')
grid on

subplot(2,2,2)
plot(c2test,h_inh(c2test,k_inh),'r-','linewidth',1.5)
hold on
plot(c2test,g_hat_test/pg_ris_ec(1),'bo','linewidth',1)
title('(b)')
xlabel('$c_2$','interpreter','latex')
ylabel('$h_2(c_2)$','interpreter','latex')
legend('true','GP')
grid on
hold off

subplot(2,2,[3,4])
plot(v,'r.','markersize',12)
title('(c)')
hold on
plot(y_hat,'bo--','linewidth',0.8)
ylim([0,1.4])
legend({'true','GP'},'location','northwest')
xlabel('$[t]$','interpreter','latex')
ylabel('$v[t]$','interpreter','latex')
grid on
hold off
