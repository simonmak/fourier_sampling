function gam_val = comp_wts(Gam_vec,w_vec,s_vec,s_max)

    % Computes Product, Order and Smoothness-Dependent (POSD) weights
    % - Gamma_vec : order wts
    % - w_vec : product wts
    % - s_vec : smoothness wts
    % - s_max : maximum smoothness

    d = length(w_vec);
    gam_mtx = permn(0:s_max,d);
    gam_val = zeros(size(gam_mtx,1),1);
    for (i = 1:size(gam_mtx,1))
        run_prod = 1.0;
        act_idx = find(gam_mtx(i,:)>0);
        for (j = act_idx)
            run_prod = run_prod * w_vec(j) * s_vec(gam_mtx(i,j)+1);
            % run_prod = run_prod * (w_vec(j)*(gam_mtx(i,j)>0.0)) * s_vec(gam_mtx(i,j)+1);
        end
        gam_val(i) = Gam_vec(nnz(gam_mtx(i,:))+1) * run_prod;
    end

end