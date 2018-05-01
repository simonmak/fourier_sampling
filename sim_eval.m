%% Simulation description:
% - Legendre / Chebyshev polynomials
% - Function evaluations (no direct Fourier sampling)

%% Simulation settings
close all
clearvars

%func_str = 'chsan10';
func_str = 'DettePepel10curv';

switch func_str
   case 'chsan10'
      % For the Cheng and Sandu Function
      d = 2;     %nominal dimension
      weights = [1 1];
      domain = [-1;1]*ones(1,d);
      func = @(xx) chsan10(xx,domain,weights);
      %func_str = 'chsan10';
      Gam_vec = ones(1,d+1); %\Gamma (order wts) for approximated function
      %Gam_vec = 1./factorial(0:d); %\Gamma (order wts) for approximated function
      s_max = 10; %maximum smoothness for approximated function
      %s_vec = 1./(( (0:s_max) + 1).^4); %s (smoothness wts) for approximated function
      s_vec = [1 1./((1:s_max).^4)]; %s (smoothness wts) for approximated function
   case 'DettePepel10curv'
      %  For the Dette
      d = 3;
      func = @(xx) detpep10curv(xx);
      Gam_vec = ones(1,d+1); %\Gamma (order wts) for approximated function
      %Gam_vec = 1./factorial(0:d); %\Gamma (order wts) for approximated function
      s_max = 7; %maximum smoothness for approximated function
      s_vec = [1 1./((1:s_max).^2)]; %s (smoothness wts) for approximated function
   otherwise
      return
end



gam_mtx = permn(0:s_max,d); %compute this just once for all
nBasis = size(gam_mtx,1); %number of basis elements
n_app = 2^14; %number of points to use for approximating L_inf norm
%basisFun = @legendreBasis; %Legendre polynomials
basisFun = @chebyshevBasis; %Chebyshev polynomials
ac_flg = true; %arc-cos flag
C = 1.5; % inflation factor
n0 = 3; % pilot sample

w_vec = -1*ones(d,1); %dummy product weights (this is not known)

% Error tolerances
num_eps = 25; %no. of errors
min_log10_eps = -3;
max_log10_eps = -1;
eps_vec = 10.^(linspace(min_log10_eps,max_log10_eps,num_eps))';
eps_vec = flipud(eps_vec); %so that we visualize the smallest tolerance

%% Set initial design and estimate Fourier coefficients
% Set desired sampling indices J for initial design
gam_mtx = permn(0:s_max,d);
samp_idx_ini = find(sum(gam_mtx ~= 0,2) <= 1 ... %no more than one nonzero element in row
  & max(gam_mtx,[],2) <= n0); %largest wavenumber is small enough
no_pts = size(samp_idx_ini,1);
[whCoord,~] = ind2sub([d,no_pts],find(gam_mtx(samp_idx_ini,:)' ~= 0));
gam_ini = gam_mtx(samp_idx_ini,:);

% Compute design and collect response
DD_ini = vdcorput_sg(gam_ini,2,ac_flg);
yy_ini = func(DD_ini);
% scatter(DD_ini(:,1),DD_ini(:,2)) %visualize in 2d

% Estimate Fourier coefficients via interpolation
[XX,basisVal_ini] = eval_X_four(DD_ini,basisFun,gam_ini,s_max);
fhat = XX\yy_ini; %least squares solution for fourier coefficients
four_coef_est_ini = zeros(nBasis,1);
four_coef_est_ini(samp_idx_ini) = fhat; %initial Fourier coefficient estimates
XX_ini = XX;
% yy - XX*fhat %check interpolation
% cond(XX'*XX) %check cond'n number

%% Construct evaluation set for testing
p = sobolset(d);
p = scramble(p,'MatousekAffineOwen');
sob_pts = 2*net(p,n_app) - 1; %stretch to fill the cube [-1,1]^d
f_true = func(sob_pts);
%four_coef = zeros(size(gam_mtx,1),1); %dummy vector -- we just want basisVal
basisValTest = eval_Basis(sob_pts,basisFun,s_max);
%test
% XX_test = eval_X_four(basisValTest,gam_mtx(samp_idx,:));

%% Run algorithm for different error tolerances
err_vec(length(eps_vec),1) = 0; %container for errors
n_vec(length(eps_vec),1) = 0; %container for sample sizes
samp_idx(nBasis,1)=0; %Want this array and those below not to change in size
DD(nBasis,d) = 0; %So we only index the elements that we need when calling them
yy(nBasis,1) = 0;
XX(nBasis,nBasis) = 0;
basisVal(nBasis,d,s_max+1) = 0;

for m = 1:length(eps_vec)
    m
    
    % Algorithm:
    % 1) Sequentially observe function to satisfy (approximate) bound:
    cur_n = no_pts; %number of points used so far
    samp_idx(1:cur_n) = samp_idx_ini; %load initial data
    four_coef_est = four_coef_est_ini; %initial Fourier coefficient estimates
    basisVal(1:cur_n,:,:) = basisVal_ini; %basis function values
    DD(1:cur_n,:) = DD_ini;
    yy(1:cur_n) = yy_ini;
    XX(1:cur_n,1:cur_n) = XX_ini;
    
    % Compute w values and sample size from pilot sample
    [n_bd,gam_val_est,w_est,~,f_hat_nm,errBd] = ...
       samp_sz(four_coef_est,Gam_vec,w_vec,s_vec,gam_mtx, ...
       eps_vec(m),C,n0,samp_idx(1:cur_n),false,false);
    f_hat_nm
    
    gam_tmp = gam_val_est; % rank remaining indices
    gam_tmp(samp_idx(1:cur_n)) = 0;  
    [elaine,rem_gam_idx] = sort(gam_tmp,'descend'); 
    rem_gam_idx = rem_gam_idx(1:end-cur_n); %remaining gamms indices to sample
    ct = 0; % counter for while loop
    
    while (errBd > eps_vec(m)) %... if sample size not enough, add largest unobserved gamma
        %disp( ['m:' num2str(m) ' - ' num2str(cur_n) '/' num2str(n_bd) ] )
        ct = ct + 1;
        %find largest unobserved gamma
        new_idx = rem_gam_idx(ct);
        cur_n = cur_n + 1;        
        samp_idx(cur_n) = new_idx;
        new_DD = vdcorput_sg(gam_mtx(new_idx,:),2,ac_flg);
        DD(cur_n,:) = new_DD;
        yy(cur_n) = func(new_DD);
        
        %recompute Fourier estimates
        [basisVal,XX] = ...
           eval_X_add(new_DD,basisVal,XX,basisFun, ...
           gam_mtx(samp_idx(1:cur_n),:),s_max,cur_n,cur_n);
        fhat = XX(1:cur_n,1:cur_n)\yy(1:cur_n); %least squares solution
        four_coef_est = zeros(nBasis,1);
        four_coef_est(samp_idx(1:cur_n)) = fhat;
        %recompute sample size bound
        [n_bd,gam_val_est,w_est,~,f_norm,errBd] = ...
           samp_sz(four_coef_est,Gam_vec,w_est,s_vec,gam_mtx, ...
            eps_vec(m),C,n0,samp_idx(1:cur_n),false,true,gam_val_est);
    end
    
    n_vec(m) = cur_n;
    f_norm

    % 2) Compute true error between f and f_app    
    %construct approximation f_app
    f_app = ...
       eval_f_four([],basisValTest,gam_mtx(samp_idx(1:cur_n),:),s_max,four_coef_est(samp_idx(1:cur_n))); 
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

rat_vec = err_vec./eps_vec; %sample size vs error ratios

% save (except large variables)
varRem = {'XX','rem_gam_idx','gam_tmp','gam_val_est','four_coef_est_ini'};
clear(varRem{:})
save(['sim_eval_results_' func_str '_ac' int2str(ac_flg) '.mat'])

sim_eval_plot_results

