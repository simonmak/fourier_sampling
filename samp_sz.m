function [nn,gam_val,w_vec,s_vec,gam_idx,f_hat_nm,errBd] = ...
   samp_sz(four_coef,Gam_vec,w_vec,w_flg,s_vec,s_flg,gam_mtx,epsTol,C,n0,nm_flg, ...
   samp_idx,gam_val)
    
   % Computes sample size given desired error tolerance eps
   % - four_coef : true or estimated fourier coefficients
   % - Gam_vec : order weights
   % - w_vec : coordinate weights
   % - w_flg : true if we know product wts
   % - s_vec : smoothness weights
   % - s_flg : true if we know smoothness wts
   % - gam_mtx : matrix of possible wavenumbers
   % - epsTol : desired error tolerance (don't use eps = machine epsilon)
   % - C : inflation factor
   % - n0 : # of pilot samples per coordinate
   % - nm_flg : true if we know ||\hat{f}||_\gamma
   % - samp_idx : index of bases in pilot sample (to compute product &
   % smoothness weights)
   % - gam_val : gamma values computed from Gam_vec, w_vec, and s_vec
   
   w_max = 1;
   sm_max = 1;
   [nBasis,d] = size(gam_mtx);
    
   % Split by flags
   if (~w_flg&&~s_flg) %Need to estimate both product and smoothness weights
      %if samp_idx not provided, sample exact coefficients for 1-d weights (up to n0 smoothness)
      if ~exist('samp_idx','var')
         samp_idx = find(sum(gam_mtx ~= 0,2) == 1 ... %just one nonzero element in row
            & max(gam_mtx,[],2) <= n0); %largest wavenumber is small enough
      else
         samp_idx(samp_idx==1) = []; %remove intercept (if in samp_idx)
      end
      nSamp_idx = size(samp_idx,1);
      [whCoord,~] = ind2sub([d,nSamp_idx],find(gam_mtx(samp_idx,:)' ~= 0));
      whSmooth = sum(gam_mtx(samp_idx,:),2);

      %Estimate product weights
      fpOvers = abs(four_coef(samp_idx)) ./ s_vec(whSmooth)';
      uell = zeros(1,d);
      for ell = 1:d
         uell(ell) = max(fpOvers(whCoord == ell));
      end
%       f_hat_pre = max(max(uell)/w_max,abs(four_coef(1))); %smallest the norm can be
      f_hat_pre = max(uell)/w_max; %smallest the norm can be
      w_vec = uell/f_hat_pre; %smallest the w_l can be
      % samp_idx = [1; samp_idx]; %put back in the intercept
      
      %Estimate smoothness weights
      fpOvers = abs(four_coef(samp_idx)) ./ w_vec(whCoord)';
      uell = zeros(1,n0);
      for ell = 1:n0
         uell(ell) = max(fpOvers(whSmooth == ell));
      end
      f_hat_pre = max(uell)/sm_max; %smallest the norm can be
      s_vec = uell/f_hat_pre; %smallest the s_j can be
      
%       abs(w_vec_tr - w_vec)
%       abs(s_vec_tr - s_vec)
      
   elseif(~w_flg&&s_flg) %Need to estimate product weights
       
      %if samp_idx not provided, sample exact coefficients for 1-d weights (up to n0 smoothness)
      if ~exist('samp_idx','var')
         samp_idx = find(sum(gam_mtx ~= 0,2) == 1 ... %just one nonzero element in row
            & max(gam_mtx,[],2) <= n0); %largest wavenumber is small enough
      else
         samp_idx(samp_idx==1) = []; %remove intercept (if in samp_idx)
      end
      nSamp_idx = size(samp_idx,1);
      [whCoord,~] = ind2sub([d,nSamp_idx],find(gam_mtx(samp_idx,:)' ~= 0));
      whSmooth = sum(gam_mtx(samp_idx,:),2);

      %Estimate product weights
      fpOvers = abs(four_coef(samp_idx)) ./ s_vec(whSmooth)';
      uell = zeros(1,d);
      for ell = 1:d
         uell(ell) = max(fpOvers(whCoord == ell));
      end
%       f_hat_pre = max(max(uell)/w_max,abs(four_coef(1))); %smallest the norm can be
      f_hat_pre = max( max(uell)/w_max,max(abs(four_coef)) ); %smallest the norm can be
      w_vec = uell/f_hat_pre; %smallest the w_l can be
      
   else
       % use provided product and smoothness weights
   end
    
   %Compute gamma weights
   if ~exist('gam_val','var')
      gam_val = comp_wts(Gam_vec,w_vec,s_vec,gam_mtx);
   end
   [gam_val_rk,gam_idx] = sort(gam_val,'descend'); %sort them in order of importance
    
   if nm_flg %if we know norm...       
      % Compute (exact) f_hat from all fourier coefficients
      f_hat_nm = max(abs(four_coef) ./ gam_val); 
   else %we do not know the norm but bound it using Fourier coefficients indexed by samp_idx   
      % Approx. norm of f_hat
      f_hat_nm = C * max( abs(four_coef([1; samp_idx])) ./ gam_val([1; samp_idx])); 
   end
    
    % Set sample size
    jout = true(nBasis,1);
    jout(1) = false;
    jout(samp_idx) = false;
    errBd = f_hat_nm * sum(gam_val(jout));
    nn = find(cumsum(gam_val_rk,'reverse') <= epsTol/f_hat_nm,1); %set sample size
    nn = max(nn,size(samp_idx,1));

end
