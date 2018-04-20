% Simulations:
% - Legendre polynomials
% - Direct Fourier sampling (no f'n evals)

%% Simulation settings
close all
clearvars
d = 5;     %dimension
Gam_vec = 1./factorial(0:d); %\Gamma (order wts)
w_vec = 1./((1:d).^2); %w (product wts) if known
s_max = 5; %maximum smoothness
s_vec = 1./(( (0:s_max) +1).^2); %s (smoothness wts)
gam_mtx = permn(0:s_max,d); %compute this just once for all
nBasis = size(gam_mtx,1) %number of basis elements
n_app = 2^14; %number of points to use for approximating L_inf norm
nm_flg = false; % do we know l-\infty norm?
w_flg = false; % do we know product weights?
rand_flg = true; % random +/- of Fourier coefficients?
randCoordOrder_flg = false; %randomize the order of the weights

if ~w_flg
    w_vec = -1*ones(1,d); %set dummy weights if no flag
end

% Error tolerances
num_eps = 25; %no. of errors
min_log10_eps = -4;
max_log10_eps = -1;
eps_vec = 10.^(linspace(min_log10_eps,max_log10_eps,num_eps))';

% Fourier coef. setting for true f'n f
Gam_vec_tr = 1./factorial(0:d); %\Gamma (order wts)
w_vec_tr = 1./((1:d).^2); %w (product wts)
s_max_tr = 5; %maximum smoothness
s_vec_tr = 1./(( (0:s_max_tr) +1).^2); %s (smoothness wts)
if randCoordOrder_flg
   w_vec_tr = w_vec_tr(randperm(d));
end
   

C = 1.2; % inflation factor
n0 = 3; % pilot sample

%% Compute true function

% Compute Fourier coef and gammas for true function
four_coef = comp_wts(Gam_vec_tr,w_vec_tr,s_vec_tr,gam_mtx);
if rand_flg
%  rad_seq = (rand(length(four_coef),1)<.5)*2 - 1; % sample Rademacher(0.5)
   rad_seq = rand(length(four_coef),1)*2 - 1; % uniform random numbers on [-1,1]
   four_coef = rad_seq .* four_coef;
end

%gam_val = comp_wts(Gam_vec,w_vec,s_vec,gam_mtx);
%[gam_val_rk,gam_idx] = sort(gam_val,'descend'); %sort for importance

p = sobolset(d);
p = scramble(p,'MatousekAffineOwen');
sob_pts = 2*net(p,n_app) - 1; %stretch to fill the cube [-1,1]^d
[f_true,basisVal,p_val] = eval_f_four(sob_pts,@legendreBasis,gam_mtx,s_max,four_coef);

%% Run algorithm for different error tolerances

err_vec(num_eps,1)=0; %container for errors
n_vec(num_eps,1)=0; %container for sample sizes

f_app(n_app,1) = 0;
for m = 1:length(eps_vec)
   %m
    
    % Algorithm:
    % 1) Compute sample size nn:
    [nn,gam_val,w_est] = samp_sz(four_coef,Gam_vec,w_vec,s_vec,gam_mtx, ...
       p_val,eps_vec(m),C,n0,[],nm_flg,w_flg);
    [gam_val_rk,gam_idx] = sort(gam_val,'descend'); 
    n_vec(m) = nn;

    % 2) Compute true error between f and f_app
    
    %evaluate f_app
    [f_app,basisVal] = ...
       eval_f_four([],basisVal,gam_mtx(gam_idx(1:nn),:),s_max,four_coef(gam_idx(1:nn)));
    
    %Record true error
    err_vec(m) = max(abs(f_true - f_app));

end

rat_vec = err_vec./eps_vec; %sample size vs error ratios
save sim_direct_results.mat %save results to plot later

%%Plot results
sim_direct_plot_results
 
