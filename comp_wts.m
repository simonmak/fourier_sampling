function gam_val = comp_wts(Gam_vec,w_vec,s_vec,gam_mtx)

    % Computes Product, Order and Smoothness-Dependent (POSD) weights
    % - Gamma_vec : order wts
    % - w_vec : product wts
    % - s_vec : smoothness wts
    % - s_max : maximum smoothness

    nBasis = size(gam_mtx,1);
    gam_val(nBasis,1) = 0;
    gam_val(1) = 1;
    for i = 2:nBasis
       act_idx = gam_mtx(i,:) > 0;
       gam_val(i) = Gam_vec(nnz(gam_mtx(i,:))+1) .* prod(w_vec(act_idx) .* s_vec(gam_mtx(i,act_idx)+1));
    end

end
