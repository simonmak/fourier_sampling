%% Simulation description:
% - Legendre / Chebyshev polynomials
% - Function evaluations (no direct Fourier sampling)

%% Simulation settings
%cd 'Fourier sampling'
close all
clearvars

%Setting the fusnction
func_str = 'chsan10';
d = 6;     %nominal dimension
n_app = 2^14; %number of points to use for approximating L_inf norm
%basisFun = @legendreBasis; %Legendre polynomials
basisFun = @chebyshevBasis; %Chebyshev polynomials
ac_flg = true; %arc-cos flag
s_flg = true; % we know smoothness weights?
%s_flg = false; % we do NOT know smoothness weights?

switch func_str
   case 'chsan10'
      plotTitle = 'Cheng \& Sandu (2010) Function';
      % For the Cheng and Sandu Function
%       weights = ones(1,d);
      weights = 1./(1:d).^3;
      weights = weights(randperm(d));
      domain = [-1;1]*ones(1,d);
      func = @(xx) chsan10(xx,domain,weights);
%       Gam_vec = ones(1,d+1); %\Gamma (order wts) for approximated function
      Gam_vec = 1./factorial(0:d); %\Gamma (order wts) for approximated function
      s_max = 5; %maximum smoothness for approximated function
      s_vec = 1./((1:s_max)); %s (smoothness wts) guess
      C = 1.1;
      
      % Error tolerances
      num_eps = 25; %no. of errors
      min_log10_eps = -3;
      max_log10_eps = -1;
      eps_vec = 10.^(linspace(min_log10_eps,max_log10_eps,num_eps))';
      eps_vec = flipud(eps_vec); %so that we visualize the smallest tolerance
   case 'DettePepel10curv'
      %  For the Dette
      func = @(xx) detpep10curv(xx);
      Gam_vec = ones(1,d+1); %\Gamma (order wts) for approximated function
      %Gam_vec = 1./factorial(0:d); %\Gamma (order wts) for approximated function
      s_max = 7; %maximum smoothness for approximated function
      s_vec = 1./((1:s_max).^2); %s (smoothness wts) guess
      
      % Error tolerances
      num_eps = 25; %no. of errors
      min_log10_eps = -3;
      max_log10_eps = -1;
      eps_vec = 10.^(linspace(min_log10_eps,max_log10_eps,num_eps))';
      eps_vec = flipud(eps_vec); %so that we visualize the smallest tolerance
    case 'otlcircuit'   
      % For the OTL circuit function
      plotTitle = 'OTL Circuit Function';
      func = @(xx) otlcircuit(xx);
%       Gam_vec = ones(1,d+1); %\Gamma (order wts) for approximated function
      Gam_vec = 1./(factorial(0:d).^2); %\Gamma (order wts) for approximated function
      s_max = 5; %maximum smoothness for approximated function
      s_vec = 1./((1:s_max)); %s (smoothness wts) guess
      C = 1.1; % inflation factor
      
      % Error tolerances
      num_eps = 25; %no. of errors
      min_log10_eps = log10(0.05);
      max_log10_eps = 0;
      eps_vec = 10.^(linspace(min_log10_eps,max_log10_eps,num_eps))';
      eps_vec = flipud(eps_vec); %so that we visualize the smallest tolerance
    
    case 'piston'   
      % For the OTL circuit function
      d = 7;
      func = @(xx) piston(xx);
%       Gam_vec = ones(1,d+1); %\Gamma (order wts) for approximated function
      Gam_vec = 1./(factorial(0:d)); %\Gamma (order wts) for approximated function
      s_max = 5; %maximum smoothness for approximated function
      s_vec = 1./((1:s_max).^4); %s (smoothness wts) guess
      C = 1.1; % inflation factor
      
      % Error tolerances
      num_eps = 25; %no. of errors
      min_log10_eps = -2;
      max_log10_eps = 0;
      eps_vec = 10.^(linspace(min_log10_eps,max_log10_eps,num_eps))';
      eps_vec = flipud(eps_vec); %so that we visualize the smallest tolerance
   otherwise
      return
end

gam_mtx = permn(0:s_max,d); %compute this just once for all
nBasis = size(gam_mtx,1); %number of basis elements
n0 = s_max; % size of pilot sample in each coordinate
w_vec = -1*ones(d,1); %dummy product weights (this is not known)

%% Set initial design and estimate Fourier coefficients
% Set desired sampling indices J for initial design
samp_idx_ini = find(sum(gam_mtx ~= 0,2) <= 1 ... %no more than one nonzero element in row
  & max(gam_mtx,[],2) <= n0); %largest wavenumber is small enough
no_pts = size(samp_idx_ini,1);
[whCoord,~] = ind2sub([d,no_pts],find(gam_mtx(samp_idx_ini,:)' ~= 0));
gam_ini = gam_mtx(samp_idx_ini,:);

% Compute design and collect response
DD_ini = vdcorput_sg(gam_ini,2,ac_flg);
yy_ini = zeros(size(gam_ini,1),1);
for (ii = 1:size(gam_ini,1))
    yy_ini(ii) = func(DD_ini(ii,:));
end
% scatter(DD_ini(:,1),DD_ini(:,2)) %visualize in 2d

% Estimate Fourier coefficients via interpolation
[XX,basisVal_ini] = eval_X_four(DD_ini,basisFun,gam_ini,s_max);
fhat = XX\yy_ini; %least squares solution for fourier coefficients
four_coef_est_ini = sparse(nBasis,1);
four_coef_est_ini(samp_idx_ini) = fhat; %initial Fourier coefficient estimates
XX_ini = XX;
% yy - XX*fhat %check interpolation
% cond(XX'*XX) %check cond'n number

%% Construct evaluation set for testing
p = sobolset(d);
p = scramble(p,'MatousekAffineOwen');
sob_pts = 2*net(p,n_app) - 1; %stretch to fill the cube [-1,1]^d
f_true = zeros(n_app,1);
for (ii = 1:n_app)
    f_true(ii) = func(sob_pts(ii,:));
end
%four_coef = zeros(size(gam_mtx,1),1); %dummy vector -- we just want basisVal
basisValTest = eval_Basis(sob_pts,basisFun,s_max);
%test
% XX_test = eval_X_four(basisValTest,gam_mtx(samp_idx,:));

%% Run algorithm for different error tolerances
err_vec(length(eps_vec),1) = 0; %container for errors
n_vec(length(eps_vec),1) = 0; %container for sample sizes
samp_idx=sparse(nBasis,1); %Want this array and those below not to change in size
DD=sparse(nBasis,d); %So we only index the elements that we need when calling them
yy=sparse(nBasis,1);
% XX(nBasis,nBasis) = 0;
XX = sparse(nBasis,nBasis);
basisVal(nBasis,d,s_max+1) = 0;

%Run algorithm
m = 1;
last_eps = eps_vec(length(eps_vec));
    
% (1) Compute w values and sample size from pilot sample
cur_n = no_pts; %number of points used so far
samp_idx(1:cur_n) = samp_idx_ini; %load initial data
four_coef_est = four_coef_est_ini; %initial Fourier coefficient estimates
basisVal(1:cur_n,:,:) = basisVal_ini; %basis function values
DD(1:cur_n,:) = DD_ini;
yy(1:cur_n) = yy_ini;
XX(1:cur_n,1:cur_n) = XX_ini;

[n_bd,gam_val_est,w_vec,s_vec,~,f_hat_nm,errBd] = ...
    samp_sz(four_coef_est,Gam_vec,w_vec,false,s_vec,s_flg, ...
    gam_mtx,last_eps,C,n0,samp_idx(1:cur_n));
f_hat_nm

gam_val_est = sparse(gam_val_est);
gam_tmp = gam_val_est; % rank remaining indices
gam_tmp(samp_idx(1:cur_n)) = 0;
[~,rem_gam_idx] = sort(gam_tmp,'descend');
rem_gam_idx = sparse(rem_gam_idx(1:end-cur_n)); %remaining gamms indices to sample
ct = 0; % counter for while loop

% (2) Sequentially observe data until adaptive error bound satisfied
while (errBd > last_eps) 
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
    %(with fixed product and smoothness weights from pilot sample)
    [n_bd,gam_val_est,w_vec,s_vec,~,f_norm,errBd] = ...
        samp_sz(four_coef_est,Gam_vec,w_vec,true,s_vec,true, ...
        gam_mtx,last_eps,C,n0,false,samp_idx(1:cur_n),gam_val_est);
    %         [n_bd,gam_val_est,w_vec,s_vec,~,f_norm,errBd] = ...
    %            samp_sz(four_coef_est,Gam_vec,w_vec,true,s_vec,false, ...
    %            gam_mtx,eps_vec(m),C,n0,samp_idx(1:cur_n),false,true,gam_val_est);
    
    % Update errors in err_vec when satisfied by adaptive error bound
    if (errBd < eps_vec(m))
        disp( ['m = ' num2str(m) ': n = ' num2str(cur_n)] )
        % Update sample size
        n_vec(m) = cur_n;
        f_norm
        % Compute true error between f and f_app    
        f_app = ...
           eval_f_four([],basisValTest,gam_mtx(samp_idx(1:cur_n),:),s_max,four_coef_est(samp_idx(1:cur_n))); 
        %Record true error
        err_vec(m) = max(abs(f_true - f_app));
        m = m + 1;
    end

end

rat_vec = err_vec./eps_vec; %sample size vs error ratios
rat_vec

% save (except large variables)
varRem = {'XX','rem_gam_idx','gam_tmp','gam_val_est','four_coef_est_ini'};
clear(varRem{:})
save(['sim_eval_results_' func_str '_d' int2str(d) '_sflg' int2str(s_flg) '.mat'])

sim_eval_plot_results(func_str,d,s_flg)

