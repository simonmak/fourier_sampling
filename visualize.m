% Assume known weights \gamma and ||\hat{f}||
% -> Legendre polynomials
% -> \hat{f}_j := \gamma_j (exact weights)

d = 2;     %dimension
s_max = 3; %maximum smoothness
Gam_vec = [1 0.5 0.5]; %\Gamma (order wts)
w_vec = [1 1]; %w (product wts)
s_vec = [1 0.5 0.5 0.5];
% s_vec = exp(-1.0*(0:s_max)); %s (smoothness wts)

eps = 0.01; %error tolerance

% compute & rank gamma
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

% Compute f_true
n_grd = 100;
x_grd = -1:(2/n_grd):1;
f_true = zeros(length(x_grd),length(x_grd));

for (i = 1:size(gam_mtx,1))
    cur_idx = gam_mtx(i,:);
    run_sum = gam_val(i)*(legendreP(cur_idx(1),x_grd))'*legendreP(cur_idx(2),x_grd);
    run_sum = run_sum / sqrt(2/(2*cur_idx(1)+1)) / sqrt(2/(2*cur_idx(2)+1)); %normalize
    f_true = f_true + run_sum;
end

%Visualize
figure
rotate3d on
[plotX,plotY] = meshgrid(x_grd);
mesh(plotX,plotY,f_true)

