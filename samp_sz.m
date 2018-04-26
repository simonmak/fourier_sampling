function [nn,gam_val,w_est,samp_idx] = samp_sz(four_coef,Gam_vec,w_vec,s_vec,gam_mtx,p_val,epsTol,C,n0,samp_idx,nm_flg,w_flg,gam_val)
    
    % Computes sample size given desired error tolerance eps
    % - four_coef : true fourier coefficients
    % - Gam_vec : order weights
    % - w_vec : product weights
    % - s_vec : smoothness weights
    % - gam_mtx : matrix of possible wavenumbers
    % - p_val : l-infty norms of basis polynomials
    % - epsTol : desired error tolerance (don't use eps = machine epsilon)
    % - C : inflation factor
    % - n0 : initial samples
    % - nm_flg : true if we know ||\hat{f}||_\gamma
    % - w_flg : true if we know product wts
    % - w_ini : initial w for optimization
   
    w_max = 1;
    [nBasis,d] = size(gam_mtx);
    
    % Set p_val to all ones if empty
    if isempty(p_val)
        p_val = ones(nBasis,1);
    end
        
    % Split by flags
    if nm_flg %do if we know norm...
        
        % Rank gamma from largest to smallest
        gam_val = comp_wts(Gam_vec,w_vec,s_vec,gam_mtx);
        gam_val_rk = sort(gam_val,'descend'); 
        
        % Compute (exact) f_hat
        f_hat_nm = max(abs(four_coef) .* p_val ./ gam_val); 
        % Set sample size
        nn = length(gam_val_rk)-sum(cumsum(flipud(gam_val_rk)) <= epsTol/f_hat_nm); 
        w_est = w_vec;
        
    elseif w_flg %do if we know weights (but not norm)...
        
        % Rank gamma from largest to smallest
        if isempty(gam_val)
            gam_val = comp_wts(Gam_vec,w_vec,s_vec,gam_mtx);
        end
        [gam_val_rk,gam_idx] = sort(gam_val,'descend'); 
        
        % Pick indices with largest n0 weights
        samp_idx = gam_idx(1:n0); 
        % Approx. norm of f_hat
        f_hat_nm = C * max( abs(four_coef(samp_idx)) .* p_val(samp_idx) ./ gam_val(samp_idx) ); 
        % Set sample size
        nn = nBasis-sum(cumsum(flipud(gam_val_rk)) <= epsTol/f_hat_nm); 
        nn = max(nn,n0);
        w_est = w_vec;
        
    else %do if we don't know weights or norm
        
        %if samp_idx not provided, sample exact coefficients for 1-d weights (up to n0 smoothness)
        if isempty(samp_idx)
           samp_idx = find(sum(gam_mtx ~= 0,2) == 1 ... %just one nonzero element in row
              & max(gam_mtx,[],2) <= n0); %largest wavenumber is small enough
           nSamp_idx = size(samp_idx,1);
           [whCoord,~] = ind2sub([d,nSamp_idx],find(gam_mtx(samp_idx,:)' ~= 0));
        else
           samp_idx(samp_idx==1) = []; %remove intercept (if in samp_idx)
           nSamp_idx = size(samp_idx,1);
           [whCoord,~] = ind2sub([d,nSamp_idx],find(gam_mtx(samp_idx,:)' ~= 0));
        end
        
        %Estimate product weights
        fpOvers = abs(four_coef(samp_idx)).* p_val(samp_idx) ...
           ./ (s_vec(gam_mtx(sub2ind([nBasis d],samp_idx, whCoord))))';
        uell = zeros(1,d);
        for ell = 1:d
           uell(ell) = max(fpOvers(whCoord == ell));
        end
%        f_hat_pre = max(uell)/w_max; %smallest the norm can be
        f_hat_pre = max(max(uell)/w_max,abs(four_coef(1)*p_val(1))); %smallest the norm can be
        w_est = uell/f_hat_pre; %smallest the w_l can be
               
        %Update gamma with estimated product weights
        gam_val = comp_wts(Gam_vec,w_est,s_vec,gam_mtx);
        gam_val_rk = sort(gam_val,'descend');

        %Compute sample size
        f_hat_nm = C * max(abs(four_coef(samp_idx)) .* p_val(samp_idx) ./ gam_val(samp_idx) ); % approx. f_hat
        nn = find(cumsum(gam_val_rk,'reverse') <= epsTol/f_hat_nm,1); %set sample size
%        nn = length(gam_val_rk)-sum(cumsum(flipud(gam_val_rk)) <= eps/f_hat_nm); %set sample size
        nn = max(nn,size(samp_idx,1));
    end

end