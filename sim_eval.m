% Simulations:
% - Legendre polynomials
% - Function evaluations (no direct Fourier sampling)

% Function specification
func = @(xx) branin(xx);

% Simulation settings
d = 2;     %dimension
Gam_vec = 1./factorial(0:d); %\Gamma (order wts)
s_max = 3; %maximum smoothness
s_vec = 1./(( (0:s_max) +1).^2); %s (smoothness wts)
ac_flg = true;

w_vec = -1*ones(d,1); %dummy weights
C = 1.2; % inflation factor
n0 = s_max; % pilot sample
w_ini = 0.25*ones(d,1); %init. for w optimization

% Error tolerances
num_eps = 25; %no. of errors
min_eps = 0.0001;
max_eps = 0.1;
eps_vec = exp(linspace(log(max_eps),log(min_eps),num_eps))';

% Set desired sampling indices J
gam_mtx = permn(0:s_max,d);
gam_idx = zeros(1,d);
for (j = 1:d)
    for (k = 1:n0)
        vc = zeros(1,d);
        vc(j) = k;
        gam_idx = [gam_idx; vc];
    end
end
[~,samp_idx] = ismember(gam_idx,gam_mtx,'rows'); %1-d sampling indices
no_pts = size(gam_idx,1);

% Compute design and collect response
DD = vdcorput_sg(gam_idx,2,ac_flg);
yy = zeros(no_pts,1);
for (i = 1:no_pts)
    yy(i) = func(DD(i,:));
end
scatter(DD(:,1),DD(:,2)) %visualize in 2d

% Interpolate and estimate fourier coefficients
XX = zeros(no_pts,no_pts);
for (i = 1:no_pts)
    for (j = 1:no_pts)
        XX(i,j) = prod(legendreP(gam_idx(j,:),DD(i,:)));
    end
end
cond(XX'*XX)
four_coef = (XX'*XX)\(XX'*yy);

% Run algorithm for different error tolerances
err_vec = zeros(length(eps_vec),1); %container for errors
n_vec = zeros(length(eps_vec),1); %container for sample sizes

for (m = 1:length(eps_vec))
    
    eps = eps_vec(m);
    
    % Algorithm:
    % 1) Compute sample size nn:
    [nn,gam_val,w_est] = samp_sz(four_coef,Gam_vec,w_vec,s_vec,s_max,[],eps,C,[],samp_idx,false,false,w_ini);
    [gam_val_rk,samp_idx] = sort(gam_val,'descend'); 
    gam_idx = gam_mtx(samp_idx(1:nn),:);
    n_vec(m) = nn;

    % 2) Compute true error between f and f_app
    
    % Load polynomials for sobol' points
    n_app = 1e4;
    if (exist(['sobol_' num2str(n_app) '_' num2str(d) '.mat']) == 2)
        %if polynomials already precomputed, then load
        load(['sobol_' num2str(n_app) '_' num2str(d) '.mat'])
    else
        %... o/w compute
        p = sobolset(d);
        p = scramble(p,'MatousekAffineOwen');
        sob_pts = (net(p,n_app)-0.5)*2;
        lp = cell(d,1);
        for (j = 1:d)
            lp{j} = zeros(n_app,s_max+1);
            for (k = 0:s_max)
                lp{j}(:,k+1) = legendreP(k,sob_pts(:,j));
            end
        end
        save(['sobol_' num2str(n_app) '_' num2str(d) '.mat'], 'sob_pts', 'lp')
    end
    
    % Evaluate true function
    f_true = zeros(n_app,1);
    for (i = 1:n_app)
        f_true(i) = func(sob_pts(i,:));
    end
    
    % Re-estimate Fourier weights, and evaluate approximation
    f_app = zeros(n_app,1);
    gam_mtx = permn(0:s_max,d);
    DD = vdcorput_sg(gam_idx,2,ac_flg);
    % scatter(DD(:,1),DD(:,2)) % visualize
    yy = zeros(nn,1);
    for (i = 1:nn)
        yy(i) = func(DD(i,:));
    end
    
    XX = zeros(nn,nn);
    for (i = 1:nn)
        for (j = 1:nn)
            XX(i,j) = prod(legendreP(gam_idx(j,:),DD(i,:)));
        end
    end
    cond(XX'*XX)
    four_coef = (XX'*XX)\(XX'*yy);
    
    for (i = 1:size(gam_mtx,1))
        cur_idx = gam_mtx(i,:);
        run_prod = ones(n_app,1);
        for (j = 1:d)
            run_prod = run_prod .* (lp{j}(:,cur_idx(j)+1));
        end
        if (ismember(i,gam_idx(1:nn)))
            f_app = f_app + run_prod * four_coef(i);
        end
    end

    %Record true error
    err_vec(m) = max(abs(f_true - f_app));

    % Observation: Estimation errors of Fourier dominates error term!
    % Try:
    % - Other functions? Smoother?
    % - Other design schemes? But how to connect to adaptive Fourier sampling?
    % - Check discrepancy between true Fourier and estimates. It should
    % hold if you get the right Fourier coefs!
    % - Maybe norm estimate not used?
end

%Plot results
rat_vec = err_vec./eps_vec; %sample size vs error ratios
cols = linspace(0,10,length(rat_vec))';
figure
scatter(rat_vec,n_vec,[],eps_vec,'x')
xlim([0 1])
xlabel({'$||f-\hat{f}||_{\infty}/\epsilon$'},'Interpreter','latex','FontSize',14)
ylabel({'Sample size $n$'},'Interpreter','latex','FontSize',14)
hcb = colorbar;
title(hcb,'\epsilon','FontSize',14)

% %Visualize (only in 2-d)
% n_grd = 100;
% x_grd = -1:(2/n_grd):1;
% f_true = zeros(length(x_grd),length(x_grd));
% f_app = zeros(length(x_grd),length(x_grd));
% 
% for (i = 1:size(gam_mtx,1))
%     cur_idx = gam_mtx(i,:);
%     run_sum = gam_val(i)*(legendreP(cur_idx(1),x_grd))'*legendreP(cur_idx(2),x_grd);
%     %run_sum = run_sum / sqrt(2/(2*cur_idx(1)+1)) / sqrt(2/(2*cur_idx(2)+1)); %normalize
%     f_true = f_true + run_sum;
%     if (ismember(i,gam_idx(1:nn)))
%         f_app = f_app + run_sum;
%     end
% end
% 
% % figure
% % rotate3d on
% % [plotX,plotY] = meshgrid(x_grd);
% % mesh(plotX,plotY,f_true)
% % 
% % figure
% % rotate3d on
% % [plotX,plotY] = meshgrid(x_grd);
% % mesh(plotX,plotY,f_app)


