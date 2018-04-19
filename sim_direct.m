% Simulations:
% - Legendre polynomials
% - Direct Fourier sampling (no f'n evals)

% Simulation settings
close all
d = 5;     %dimension
Gam_vec = 1./factorial(0:d); %\Gamma (order wts)
w_vec = 1./((1:d).^2); %w (product wts) if known
s_max = 3; %maximum smoothness
s_vec = 1./(( (0:s_max) +1).^2); %s (smoothness wts)
gam_mtx = permn(0:s_max,d); %compute this just once for all
nm_flg = false; % do we know l-\infty norm?
w_flg = false; % do we know product weights?
rand_flg = true; % random +/- of Fourier coefficients?
randCoordOrder_flg = false; %randomize the order of the weights

if ~w_flg
    w_vec = -1*ones(1,d); %set dummy weights if no flag
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
if randCoordOrder_flg
   w_vec_tr = w_vec_tr(randperm(d));
end
   

C = 1.2; % inflation factor
n0 = s_max; % pilot sample

%w_ini = 0.25*ones(1,d); %init. for w optimization

% Compute Fourier coef and gammas
four_coef = comp_wts(Gam_vec_tr,w_vec_tr,s_vec_tr,gam_mtx);
if rand_flg
%  rad_seq = (rand(length(four_coef),1)<.5)*2 - 1; % sample Rademacher(0.5)
   rad_seq = rand(length(four_coef),1)*2 - 1; % uniform random numbers on [-1,1]
   four_coef = rad_seq .* four_coef;
end

gam_val = comp_wts(Gam_vec,w_vec,s_vec,gam_mtx);
[gam_val_rk,gam_idx] = sort(gam_val,'descend'); %sort for importance

% Run algorithm for different error tolerances

err_vec = zeros(length(eps_vec),1); %container for errors
n_vec = zeros(length(eps_vec),1); %container for sample sizes

for m = 1:length(eps_vec)
    
    eps = eps_vec(m);
    
    % Algorithm:
    % 1) Compute sample size nn:
    [nn,gam_val,w_est] = samp_sz(four_coef,Gam_vec,w_vec,s_vec,gam_mtx,[],eps,C,n0,[],nm_flg,w_flg,w_ini);
    [gam_val_rk,gam_idx] = sort(gam_val,'descend'); 
    n_vec(m) = nn;

    % 2) Compute true error between f and f_app
    n_app = 1e4;
    if exist(['sobol_' num2str(n_app) '_' num2str(d) '.mat'],'file')
        %if polynomials already precomputed, then load
        load(['sobol_' num2str(n_app) '_' num2str(d) '.mat'])
    else
        %... o/w compute
        p = sobolset(d);
        p = scramble(p,'MatousekAffineOwen');
        sob_pts = net(p,n_app);
        lp = cell(d,1);
        for j = 1:d
            lp{j} = zeros(n_app,s_max+1);
            for k = 0:s_max
                lp{j}(:,k+1) = legendreP(k,sob_pts(:,j));
            end
        end
        save(['sobol_' num2str(n_app) '_' num2str(d) '.mat'], 'sob_pts', 'lp')
    end

    f_true = zeros(n_app,1);
    f_app = zeros(n_app,1);
    gam_mtx = permn(0:s_max,d);
    
    %evaluate f and f_app
    for i = 1:size(gam_mtx,1)
        cur_idx = gam_mtx(i,:);
        run_prod = ones(n_app,1);
        for j = 1:d
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

rat_vec = err_vec./eps_vec; %sample size vs error ratios
save sim_direct_results.mat %save results to plot later

%%Plot results
sim_direct_plot_results
 
