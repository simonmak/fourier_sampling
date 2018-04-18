function DD = vdcorput_sg(gam_idx,b,acflg)
    % Computes (very sparse) van der Corput sparse grid from indices
    % gam_idx
    % - gam_idx : Indices to estimate
    % - b : base for vdC sequence
    % - acflg : arc-cos flag
    
    no_idx = size(gam_idx,1);
    mx = max(max(gam_idx));
    vdc_seq = vdcorput(mx,b);
    
    DD = zeros(no_idx,size(gam_idx,2));
    if (~acflg) %no arc-cos adjustment
        for (i = 1:no_idx)
            DD(i,:) = 2*mod( vdc_seq(gam_idx(i,:)+1) + 1/3, 1 )-1;
        end
    else %arc-cos adjustment
        for (i = 1:no_idx)
            DD(i,:) = -cos(pi*mod( vdc_seq(gam_idx(i,:)+1) + 1/3 , 1 ));
        end        
    end
    
end