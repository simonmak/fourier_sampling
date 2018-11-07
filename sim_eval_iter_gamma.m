%% Simulation description:
% - Legendre / Chebyshev polynomials
% - Function evaluations (no direct Fourier sampling)

%% Simulation settings
close all
clearvars

%func_str = 'chsan10';
%func_str = 'DettePepel10curv';
func_str = 'otlcircuit';

switch func_str
   case 'chsan10'
      % For the Cheng and Sandu Function
      d = 6;     %nominal dimension
      whVar = [false true true false false false];
      domain = [-1;1]*ones(1,d);
      func = @(xx) chsan10(xx,domain,whVar);
      %func_str = 'chsan10';
      Gam_vec = ones(1,d+1); %\Gamma (order wts) for approximated function
      %Gam_vec = 1./factorial(0:d); %\Gamma (order wts) for approximated function
      s_max = 30; %maximum smoothness for approximated function
      %s_vec = 1./(( (0:s_max) + 1).^4); %s (smoothness wts) for approximated function
      s_power = 4;
      num_eps = 25; %no. of errors
      min_log10_eps = -3;
      max_log10_eps = -1;

   case 'DettePepel10curv'
      %  For the Dette
      d = 3;
      func = @(xx) detpep10curv(xx);
      Gam_vec = ones(1,d+1); %\Gamma (order wts) for approximated function
      %Gam_vec = 1./factorial(0:d); %\Gamma (order wts) for approximated function
      s_max = 20; %maximum smoothness for approximated function
      s_power = 2;
      num_eps = 10; %no. of errors
      min_log10_eps = -2;
      max_log10_eps = -1;

   case 'otlcircuit'
      d = 6;
      func = @(xx) otlcircuit(xx);
      Gam_vec = ones(1,d+1); %\Gamma (order wts) for approximated function
      %Gam_vec = 1./factorial(0:d); %\Gamma (order wts) for approximated function
      s_max = 30; %maximum smoothness for approximated function
      s_power = 6;
      num_eps = 25; %no. of errors
      min_log10_eps = -1.7;
      max_log10_eps = -0;
      
   otherwise
      return
end

maxCond = 1;
s_vec = [1 1./((1:s_max).^s_power)]; %s (smoothness wts) for approximated function
s_sum = zeta(s_power);

n_app = 2^14; %number of points to use for approximating L_inf norm
%basisFun = @legendreBasis; %Legendre polynomials
basisFun = @chebyshevBasis; %Chebyshev polynomials
ac_flg = true; %arc-cos flag
C = 1.5; % inflation factor
n0 = 3; % pilot sample

w_vec = -1*ones(d,1); %dummy product weights (this is not known)

% Error tolerances
eps_vec = 10.^(linspace(min_log10_eps,max_log10_eps,num_eps))';
eps_vec = flipud(eps_vec); %so that we visualize the smallest tolerance

%% Set initial design and estimate Fourier coefficients
% Set desired sampling indices J for initial design
nBasis = d*n0+1;
nBLength = 1e3;
gam_mtx = zeros(nBLength,d);
whCoord = zeros(nBLength,1);
for k = 1:d
   gam_mtx((k-1)*n0+2:k*n0+1,k) = (1:n0)';
   whCoord((k-1)*n0+2:k*n0+1) = k;
end
no_pts = nBasis;

% Compute design and collect response
DD(nBLength,d) = 0; %So we only index the elements that we need when calling them
yy(nBLength,1) = 0;
XX(nBLength,nBLength) = 0;
basisVal(nBLength,d,s_max+1) = 0;


DD(1:no_pts,:) = vdcorput_sg(gam_mtx(1:no_pts,:),2,ac_flg);
yy(1:no_pts) = func(DD(1:no_pts,:));

% Estimate Fourier coefficients via interpolation
[XX(1:no_pts,1:no_pts),basisVal(1:no_pts,:,:)] = ...
   eval_X_four(DD(1:no_pts,:),basisFun,gam_mtx(1:no_pts,:),s_max);
fhat = XX(1:no_pts,1:no_pts)\yy(1:no_pts); %least squares solution for fourier coefficients
four_coef_est = zeros(nBLength,1);
four_coef_est(1:no_pts) = fhat; %initial Fourier coefficient estimates
maxCond = max(maxCond,cond(XX(1:no_pts,1:no_pts)));
% yy - XX*fhat %check interpolation
% cond(XX'*XX) %check cond'n number

%% Construct evaluation set for testing
p = sobolset(d);
p = scramble(p,'MatousekAffineOwen');
sob_pts = 2*net(p,n_app) - 1; %stretch to fill the cube [-1,1]^d
f_true = func(sob_pts);
basisValTest = eval_Basis(sob_pts,basisFun,s_max);

%% Run algorithm for different error tolerances
err_vec(length(eps_vec),1) = 0; %container for errors
n_vec(length(eps_vec),1) = 0; %container for sample sizes
samp_idx(nBasis,1)=0; %Want this array and those below not to change in size
gam_val(nBLength,1) = 0;
gam_idx(nBLength,1) = 0;
gam_bigs(nBLength,1) = false;
wh_w_inv(1,d)= 0;

for m = 1:length(eps_vec)
    m
    
    % Algorithm:
    % 1) Sequentially observe function to satisfy (approximate) bound:
    cur_n = no_pts; %number of points used so far
    
    % Compute w values and sample size from pilot sample
    [gam_val(1:cur_n),w_est,gam_idx(1:cur_n),f_hat_nm,err_Bd] = ...
       samp_sz_iter(four_coef_est(1:cur_n),Gam_vec,w_vec,s_power,s_sum, ...
       gam_mtx(1:cur_n,:),C,n0,(1:cur_n)',false,false);
%     f_hat_nm, err_Bd
    
    [~,wh_w_est] = sort(w_est,'descend'); %sort the coordinate weights
    wh_w_inv(wh_w_est) = 1:d;
    gamma_sum = prod(1 + w_est*s_sum);
    cur_n_gam = cur_n;
    for k = 1:d
       %add one with s incremented by one
       cur_n_gam = cur_n_gam + 1; %one addition to the wavenumber pool
       gam_mtx(cur_n_gam,:) = zeros(1,d); %add new wavenumber
       gam_mtx(cur_n_gam,k) = n0 + 1;
       gam_val(cur_n_gam) = ... %compute gamma value for this new wavenumber
          comp_wts_iter(Gam_vec,w_est,s_power,gam_mtx(cur_n_gam,:));
       %gam_idx points to where the gammas are in order of largest to
       %smallest starting with those after cur_n
       wh_insert = find(gam_val(cur_n_gam) > gam_val(gam_idx(cur_n+1:cur_n_gam-1)),1,'first'); %find where this fits in the order
       if ~isempty(wh_insert)
          gam_idx(cur_n+wh_insert:cur_n_gam) = [cur_n_gam; gam_idx(cur_n+wh_insert:cur_n_gam-1)];
       else
          gam_idx(cur_n_gam) = cur_n_gam;
       end
       for kk = wh_w_inv(k)+1:d
          for j = 1:n0
             %add additional interaction
             cur_n_gam = cur_n_gam + 1; %one addition to the wavenumber pool
             gam_mtx(cur_n_gam,:) = zeros(1,d);
             gam_mtx(cur_n_gam,k) = j;
             gam_mtx(cur_n_gam,wh_w_est(kk)) = 1; %add lowest order interaction
             gam_val(cur_n_gam) = ...
                comp_wts_iter(Gam_vec,w_est,s_power,gam_mtx(cur_n_gam,:));
             wh_insert = find(gam_val(cur_n_gam) > gam_val(gam_idx(cur_n+1:cur_n_gam-1)),1,'first'); %find where this fits in the order
             if ~isempty(wh_insert)
                gam_idx(cur_n+wh_insert:cur_n_gam) = [cur_n_gam; gam_idx(cur_n+wh_insert:cur_n_gam-1)];
             else
                gam_idx(cur_n_gam) = cur_n_gam;
             end
          end
       end
    end
    
    while (err_Bd > eps_vec(m)) %... if sample size not enough, add largest unobserved gamma
        %disp( ['m:' num2str(m) ' - ' num2str(cur_n) '/' num2str(n_bd) ] )
        cur_n = cur_n + 1; %next sample
        new_idx = gam_idx(cur_n);
        new_wavenum = gam_mtx(new_idx,:);
        new_DD = vdcorput_sg(new_wavenum,2,ac_flg);
        DD(cur_n,:) = new_DD;
        yy(cur_n) = func(new_DD);
        
        %recompute Fourier estimates and update error bound
        [basisVal,XX] = ...
           eval_X_add(new_DD,basisVal,XX,basisFun, ...
           gam_mtx(gam_idx(1:cur_n),:),s_max,cur_n,cur_n);
        fhat = XX(1:cur_n,1:cur_n)\yy(1:cur_n); %least squares solution
        four_coef_est(1:cur_n) = fhat;
        maxCond = max(maxCond,cond(XX(1:cur_n,1:cur_n)));
        %recompute error bound
        wh_nz = find([0; four_coef_est(2:cur_n) ~= 0]);
        f_hat_nm = C * max( abs(four_coef_est(wh_nz)) ...
           ./ gam_val(gam_idx(wh_nz))); 
        tail_sum = gamma_sum - sum(gam_val(gam_idx(1:cur_n)));
        err_Bd = f_hat_nm * tail_sum;
                 
        %add next two gammas to the pool
        [cur_n_gam,gam_mtx,gam_val,gam_idx] = ...
           add_new_wavenum(new_wavenum,cur_n,cur_n_gam,gam_mtx,gam_val,gam_idx, ...
              wh_w_est,wh_w_inv,Gam_vec,w_est,s_power);

%        f_hat_nm, err_Bd, tail_sum
    end
    
    n_vec(m) = cur_n;
    %f_hat_nm, cur_n

    % 2) Compute true error between f and f_app    
    %construct approximation f_app
    f_app = ...
       eval_f_four([],basisValTest,gam_mtx(gam_idx(1:cur_n),:),s_max,four_coef_est(1:cur_n)); 
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
maxCond

rat_vec = err_vec./eps_vec; %sample size vs error ratios

% save (except large variables)
varRem = {'XX','rem_gam_idx','gam_tmp','gam_val_est','four_coef_est_ini'};
clear(varRem{:})
save(['sim_eval_results_' func_str '_ac' int2str(ac_flg) '.mat'])

sim_eval_iter_plot_results(func_str,ac_flg)

