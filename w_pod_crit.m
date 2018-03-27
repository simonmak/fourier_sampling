function ret = w_pod_crit(ww, idx, Gam_vec, s_vec, s_max, four_coef)
    % Optimization criterion for estimating product weights

    gam_val = comp_wts(Gam_vec,ww,s_vec,s_max);

%     ret = max( abs(four_coef(idx)) ./ gam_val(idx) );
    ret = sum( ( abs(four_coef(idx)) ./ gam_val(idx) - 1 ).^2 );
end