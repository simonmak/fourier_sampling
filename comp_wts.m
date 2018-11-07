function gam_val = comp_wts(Gam_vec,w_vec,s_vec,gam_mtx)

    % Computes Product, Order and Smoothness-Dependent (POSD) weights
    % - Gamma_vec : order wts
    % - w_vec : product wts
    % - s_vec : smoothness wts
    % - s_max : maximum smoothness

    nBasis = size(gam_mtx,1);
    gam_val = ones(nBasis,1);
    for i = 1:nBasis
       act_idx = gam_mtx(i,:) > 0;
       if sum(act_idx)
          gam_val(i) = Gam_vec(nnz(gam_mtx(i,:))+1) .* prod(w_vec(act_idx) ...
             .* s_vec(gam_mtx(i,act_idx)));
       end
    end

end
