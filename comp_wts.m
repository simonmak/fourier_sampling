function gam_val = comp_wts(Gam_vec,w_vec,s_vec,s_max)

    % Computes Product, Order and Smoothness-Dependent (POSD) weights
    % - Gamma_vec : order wts
    % - w_vec : product wts
    % - s_vec : smoothness wts
    % - s_max : maximum smoothness

    d = length(w_vec); %dimension of input
    gam_mtx = permn(0:s_max,d);
    nbasis = size(gam_mtx,1);
    gam_val = ones(nbasis,1);
    for i = 2:nbasis
       act_idx = gam_mtx(i,:) > 0;
       gam_val(i) = Gam_vec(nnz(gam_mtx(i,:))+1) .* prod(w_vec(act_idx) .* s_vec(gam_mtx(i,act_idx)+1));
    end

end
