%%Plot results
clearvars

load sim_direct_results.mat %load results for plotting
InitializeDisplay %add some variables for nice plotting

figure
log10epsVec = log10(eps_vec);
scatter(rat_vec,n_vec,800,log10epsVec,'.')
set(gca,'YScale','log')
xlim([0 max([rat_vec*1.2;1])])
ylim(10.^[floor(log10(min(n_vec))) ceil(log10(max(n_vec)))])
xlabel({'\(||f-\hat{f}||_{\infty}/\epsilon\)'})
ylabel({'Sample size \(n\)'})
hcb = colorbar;
title(hcb,'\(\epsilon\)','Interpreter','latex')
tickVals = floor(min_log10_eps:max_log10_eps);
tickLabels = 10.^tickVals;
set(hcb,'Ticks',tickVals,'TickLabels',tickLabels, ...
   'Limits',[tickVals(1) tickVals(end)])

%% Visualize (only a 2-d projection)
xcoord = 1;
ycoord = 2;
nom_Val = 0;
n_grd = 100;
x_grd = -1:(2/n_grd):1;
n_grd_val = length(x_grd);
[xx,yy] = meshgrid(x_grd);
n_Vis = n_grd_val.^2;
x_Vis = nom_Val*ones(n_Vis,d);
x_Vis(:,xcoord) = xx(:);
x_Vis(:,ycoord) = yy(:);
[f_true_Vis,basisVal] = eval_f_four(x_Vis,@legendreBasis,gam_mtx,s_max,four_coef);
[f_app_Vis] = ...
       eval_f_four([],basisVal,gam_mtx(gam_idx(1:nn),:),s_max,four_coef(gam_idx(1:nn)));

figure
rotate3d on
f_true_Vis = reshape(f_true_Vis,[n_grd_val n_grd_val]);
surf(x_grd,x_grd,f_true_Vis); shading interp

figure
rotate3d on
f_app_Vis = reshape(f_app_Vis,[n_grd_val n_grd_val]);
surf(x_grd,x_grd,f_app_Vis); shading interp

figure(1)

