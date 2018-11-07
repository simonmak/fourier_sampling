% Simulations:
% - Legendre polynomials
% - Direct Fourier sampling (no f'n evals)

%% Simulation settings
close all
clearvars

d = 8;     %dimension

Gam_vec_tr = 1./factorial(0:d); %\Gamma (order wts) for true function
Gam_vec = 1./factorial(0:d); %\Gamma (order wts) guess

w_vec_tr = 1./((1:d).^2); %w (product wts) for true function
w_vec = -1*ones(1,d); %w (product wts) guess

s_max = 5; %maximum smoothness for approximated function
s_max_tr = s_max; %maximum smoothness for true function, I think that this has to be the same
s_vec_tr = 1./((1:s_max_tr).^3); %s (smoothness wts) for true function
% s_vec = s_vec_tr;
s_vec = 1./(1:s_max); %s (smoothness wts) guess

gam_mtx = permn(0:s_max,d); %compute this just once for all
nBasis = size(gam_mtx,1); %number of basis elements
n_app = 2^14; %number of points to use for approximating L_inf norm
basisFun = @legendreBasis; %Legendre polynomials
%basisFun = @chebyshevBasis; %Chebyshev polynomials

nm_flg = false; % do we know l-\infty norm?
w_flg = false; % do we know product weights?
s_flg = false; % do we know smoothness weights?
rand_flg = true; % random +/- of Fourier coefficients?

randCoordOrder_flg = true; %randomize the order of the weights
if randCoordOrder_flg
   w_vec_tr = w_vec_tr(randperm(d));
end

% Error tolerances
num_eps = 25; %no. of errors
min_log10_eps = -4;
max_log10_eps = -1;
eps_vec = 10.^(linspace(min_log10_eps,max_log10_eps,num_eps))';
eps_vec = flipud(eps_vec); %so that we visualize the smallest tolerance

C = 1.2; % inflation factor
n0 = s_max; % pilot sample in each coordinate

%% Compute true function

% Compute Fourier coef and gammas for true function
four_coef = comp_wts(Gam_vec_tr,w_vec_tr,s_vec_tr,gam_mtx);
if rand_flg
%  rad_seq = (rand(length(four_coef),1)<.5)*2 - 1; % sample Rademacher(0.5)
   rad_seq = rand(length(four_coef),1)*2 - 1; % uniform random numbers on [-1,1]
   four_coef = rad_seq .* four_coef;
end

p = sobolset(d);
p = scramble(p,'MatousekAffineOwen');
sob_pts = 2*net(p,n_app) - 1; %stretch to fill the cube [-1,1]^d
[f_true,basisVal] = eval_f_four(sob_pts,basisFun,gam_mtx,s_max_tr,four_coef);

%% Run algorithm for different error tolerances

err_vec(num_eps,1)=0; %container for errors
four_ignored(num_eps,1)=0; %container for L^1 norm of Fourier coefficients not used
normf_gamma_ignored(num_eps,1)=0; %container for upper bound on error
n_vec(num_eps,1)=0; %container for sample sizes

f_app(n_app,1) = 0;
for m = 1:length(eps_vec)
    m
    
    % Algorithm:
    % 1) Compute sample size nn:
    [nn,gam_val,w_vec,s_vec,gam_idx,f_hat_nm] = ...
        samp_sz(four_coef,Gam_vec,w_vec,w_flg,s_vec,s_flg,gam_mtx,eps_vec(m),C,n0,nm_flg);
    wh_gam_in = gam_idx(1:nn);
    wh_gam_out = gam_idx(nn+1:nBasis);
    n_vec(m) = nn;

    % 2) Compute true error between f and f_app    
    %evaluate f_app
    [f_app,basisVal] = ...
       eval_f_four([],basisVal,gam_mtx(wh_gam_in,:),s_max,four_coef(wh_gam_in)); 
    %Record true error
    err_vec(m) = max(abs(f_true - f_app));
    four_ignored(m) = sum(abs(four_coef(wh_gam_out)));
    normf_gamma_ignored(m) = f_hat_nm*sum(abs(gam_val(wh_gam_out)));

end

rat_vec = err_vec./eps_vec; %sample size vs error ratios
rat_four_eps_vec = four_ignored./eps_vec; %sample size vs error ratios
rat_normf_eps_vec = normf_gamma_ignored./eps_vec; %sample size vs error ratios
% save sim_direct_results.mat %save results to plot later
save(['sim_eval_results_d' int2str(d) '_sflg' int2str(s_flg) '.mat'])

%%Plot results
sim_direct_plot_results(d,s_flg)
 
