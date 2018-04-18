function [nn,gam_val,w_est] = samp_sz(four_coef,Gam_vec,w_vec,s_vec,s_max,p_val,eps,C,n0,samp_idx,nm_flg,w_flg,w_ini)
    
    % Computes sample size given desired error tolerance eps
    % - four_coef : true fourier coefficients
    % - Gam_vec : order weights
    % - w_vec : product weights
    % - s_vec : smoothness weights
    % - s_max : max. smoothness
    % - p_val : l-infty norms of basis polynomials
    % - eps : desired error tolerance
    % - C : inflation factor
    % - n0 : initial samples
    % - nm_flg : true if we know ||\hat{f}||_\gamma
    % - w_flg : true if we know product wts
    % - w_ini : initial w for optimization
   
    d = length(w_vec);
    gam_mtx = permn(0:s_max,d);
    nbasis = size(gam_mtx,1);
    
    % Set p_val to all ones if empty
    if isempty(p_val)
        p_val = ones(size(gam_mtx,1),1);
    end
        
    % Split by flags
    if nm_flg %do if we know norm...
        
        % Rank gamma from largest to smallest
        gam_val = comp_wts(Gam_vec,w_vec,s_vec,s_max);
        gam_val_rk = sort(gam_val,'descend'); 
        
        % Compute (exact) f_hat
        f_hat_nm = max(abs(four_coef) .* p_val ./ gam_val); 
        % Set sample size
        nn = length(gam_val_rk)-sum(cumsum(flipud(gam_val_rk)) <= eps/f_hat_nm); 
        w_est = w_vec;
        
    elseif w_flg %do if we know weights (but not norm)...
        
        % Rank gamma from largest to smallest
        gam_val = comp_wts(Gam_vec,w_vec,s_vec,s_max);
        [gam_val_rk,gam_idx] = sort(gam_val,'descend'); 
        
        % Pick indices with largest n0 weights
        samp_idx = gam_idx(1:n0); 
        % Approx. norm of f_hat
        f_hat_nm = C * max( abs(four_coef(samp_idx)) .* p_val(samp_idx) ./ gam_val(samp_idx) ); 
        % Set sample size
        nn = nbasis-sum(cumsum(flipud(gam_val_rk)) <= eps/f_hat_nm); 
        nn = max(nn,n0);
        w_est = w_vec;
        
    else %do if we don't know weights or norm
        
        %if samp_idx not provided, sample exact coefficients for 1-d weights (up to n0 smoothness)
        if isempty(samp_idx)
            samps = zeros(n0,d); %container for samples
            oned = zeros(1,d);
            for j = 1:d
                for k = 1:n0
                    vc = zeros(1,d);
                    vc(j) = k;
                    oned = [oned; vc];
                end
            end
            [~,samp_idx] = ismember(oned,gam_mtx,'rows'); %index for 1-d weights
            %Estimate product weights via nonlinear optimization
            fun = @(ww) w_pod_crit(ww, samp_idx, Gam_vec, s_vec, s_max, four_coef(samp_idx));
            w_est = fmincon(fun,w_ini,[],[],[],[],zeros(d,1),ones(d,1));
            
            %Update gamma with estimated product weights
            gam_val = comp_wts(Gam_vec,w_est,s_vec,s_max);
            [gam_val_rk,gam_idx] = sort(gam_val,'descend');

            %Compute sample size
            f_hat_nm = C * max( abs(four_coef(samp_idx)) .* p_val(samp_idx) ./ gam_val(samp_idx) ); % approx. f_hat
            nn = length(gam_val_rk)-sum(cumsum(flipud(gam_val_rk)) <= eps/f_hat_nm); %set sample size
            nn = max(nn,n0);
        
        else
            %Estimate product weights via nonlinear optimization
            fun = @(ww) w_pod_crit(ww, samp_idx, Gam_vec, s_vec, s_max, four_coef);
            w_est = fmincon(fun,w_ini,[],[],[],[],zeros(d,1),ones(d,1));
            
            %Update gamma with estimated product weights
            gam_val = comp_wts(Gam_vec,w_est,s_vec,s_max);
            [gam_val_rk,gam_idx] = sort(gam_val,'descend');

            %Compute sample size
            f_hat_nm = C * max( abs(four_coef) .* p_val(samp_idx) ./ gam_val(samp_idx) ); % approx. f_hat
            nn = length(gam_val_rk)-sum(cumsum(flipud(gam_val_rk)) <= eps/f_hat_nm); %set sample size
            nn = max(nn,size(samp_idx,1));
        end
        
    end

end