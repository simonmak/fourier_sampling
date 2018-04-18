% Simulations:
% - Legendre polynomials
% - Direct Fourier sampling (no f'n evals)

% Simulation settings
d = 5;     %dimension
Gam_vec = 1./factorial(0:d); %\Gamma (order wts)
w_vec = 1./((1:d).^2); %w (product wts)
s_max = 3; %maximum smoothness
s_vec = 1./(( (0:s_max) +1).^2); %s (smoothness wts)
nm_flg = false; % do we know l-\infty norm?
w_flg = false; % do we know product weights?
rand_flg = false; % random +/- of Fourier coefficients?

if (~w_flg)
    w_vec = -1*ones(d,1); %set dummy weights if no flag
end

% Error tolerances
num_eps = 25; %no. of errors
min_eps = 0.0001;
max_eps = 0.1;
eps_vec = exp(linspace(log(max_eps),log(min_eps),num_eps))';

% Fourier coef. setting for true f'n f
Gam_vec_tr = 1./factorial(0:d); %\Gamma (order wts)
w_vec_tr = 1./((1:d).^2); %w (product wts)
s_max_tr = 3; %maximum smoothness
s_vec_tr = 1./(( (0:s_max_tr) +1).^2); %s (smoothness wts)

C = 1.2; % inflation factor
n0 = s_max; % pilot sample

w_ini = 0.25*ones(d,1); %init. for w optimization

% Compute Fourier coef and gammas
four_coef = comp_wts(Gam_vec_tr,w_vec_tr,s_vec_tr,s_max_tr);
if (rand_flg)
    rad_seq = (rand(length(four_coef),1)<.5)*2 - 1; % sample Rademacher(0.5)
else
    rad_seq = ones(length(four_coef),1);
end
four_coef = rad_seq .* four_coef;
gam_val = comp_wts(Gam_vec,w_vec,s_vec,s_max);
[gam_val_rk,gam_idx] = sort(gam_val,'descend'); %sort for importance

% Run algorithm for different error tolerances

err_vec = zeros(length(eps_vec),1); %container for errors
n_vec = zeros(length(eps_vec),1); %container for sample sizes

for (m = 1:length(eps_vec))
    
    eps = eps_vec(m);
    
    % Algorithm:
    % 1) Compute sample size nn:
    [nn,gam_val,w_est] = samp_sz(four_coef,Gam_vec,w_vec,s_vec,s_max,[],eps,C,n0,[],nm_flg,w_flg,w_ini);
    [gam_val_rk,gam_idx] = sort(gam_val,'descend'); 
    n_vec(m) = nn;

    % 2) Compute true error between f and f_app
    n_app = 1e4;
    if (exist(['sobol_' num2str(n_app) '_' num2str(d) '.mat']) == 2)
        %if polynomials already precomputed, then load
        load(['sobol_' num2str(n_app) '_' num2str(d) '.mat'])
    else
        %... o/w compute
        p = sobolset(d);
        p = scramble(p,'MatousekAffineOwen');
        sob_pts = net(p,n_app);
        lp = cell(d,1);
        for (j = 1:d)
            lp{j} = zeros(n_app,s_max+1);
            for (k = 0:s_max)
                lp{j}(:,k+1) = legendreP(k,sob_pts(:,j));
            end
        end
        save(['sobol_' num2str(n_app) '_' num2str(d) '.mat'], 'sob_pts', 'lp')
    end

    f_true = zeros(n_app,1);
    f_app = zeros(n_app,1);
    gam_mtx = permn(0:s_max,d);
    
    %evaluate f and f_app
    for (i = 1:size(gam_mtx,1))
        cur_idx = gam_mtx(i,:);
        run_prod = ones(n_app,1);
        for (j = 1:d)
            run_prod = run_prod .* (lp{j}(:,cur_idx(j)+1));
        end
        f_true = f_true + run_prod * four_coef(i);
        if (ismember(i,gam_idx(1:nn)))
            f_app = f_app + run_prod * four_coef(i);
        end
    end

    %Record true error
    err_vec(m) = max(abs(f_true - f_app));

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


