%%Plot results
function sim_eval_plot_results(func_str,ac_flg)

%func_str = 'chsan10'; %function string
%ac_flg = 1; %arc-cos flag
load(['sim_eval_results_' func_str '_ac' int2str(ac_flg) '.mat'], ...
   'rat_vec', 'eps_vec', 'n_vec', 'w_est', 'func', 'd', 'cur_n', ...
   'basisFun', 'gam_mtx', 's_max', 'four_coef_est', 'samp_idx', ...
   'min_log10_eps', 'max_log10_eps') %load results for plotting
InitializeDisplay %add some variables for nice plotting

figure
log10epsVec = log10(eps_vec);
h = scatter(rat_vec,n_vec,800,log10epsVec,'.'); %plot ratio of actual error to tolerance, with color corresonding to tolerance
set(gca,'XScale','log', 'YScale','log')
xlim([min([[rat_vec ]*0.8; 0.1]) ...
 max([[rat_vec ]*1.2; 1])])
ylim(10.^[floor(log10(min(n_vec))) ceil(log10(max(n_vec)))])
%xlabel({'\(||f-\hat{f}||_{\infty}/\varepsilon\)'})
ylabel({'Sample size \(n\)'})
hcb = colorbar; %showing tolerance values
title(hcb,'\(\varepsilon\)','interpreter','latex')
tickVals = floor(min_log10_eps:max_log10_eps);
tickLabels = 10.^tickVals;
set(hcb,'Ticks',tickVals,'TickLabels',tickLabels, ...
   'Limits',[tickVals(1) tickVals(end)])
[~,leg_icons] = legend(h,{'\(||f-f_{\mbox{app}}||_{\infty}/\varepsilon\) \quad'}, ...
   'box','off','location','south','orientation','horizontal');
title( [func_str ' \((d = ' num2str(d) ')\)'] )
set(gcf,'Position',[200,200,1000,500]) %make figure big enough and the right aspect ratio
leg_icons(2).Children.MarkerSize = 30; %make legend icons large enough

%% Visualize (only a 2-d projection)
[~,whCoordBig] = sort(w_est,'descend');
xcoord = whCoordBig(1); % coordinate with the biggest w
ycoord = whCoordBig(2); % coordinate with the next biggest w
nom_Val = 0;
n_grd = 100;
x_grd = -1:(2/n_grd):1;
n_grd_val = length(x_grd);
[xx,yy] = meshgrid(x_grd);
n_Vis = n_grd_val.^2;
x_Vis = nom_Val*ones(n_Vis,d);
x_Vis(:,xcoord) = xx(:);
x_Vis(:,ycoord) = yy(:);
f_true_Vis = func(x_Vis);
[f_app_Vis] = ...
       eval_f_four(x_Vis,basisFun,gam_mtx(samp_idx(1:cur_n),:), ...
       s_max,four_coef_est(samp_idx(1:cur_n)));
zl = [min([f_app_Vis; f_true_Vis]) max([f_app_Vis; f_true_Vis]) ]; %zlims
   
figure
rotate3d on
f_true_Vis = reshape(f_true_Vis,[n_grd_val n_grd_val]);
surf(x_grd,x_grd,f_true_Vis); shading interp
zlim(zl)
xlabel(['\(x_{' int2str(whCoordBig(1)) '}\)'])
ylabel(['\(x_{' int2str(whCoordBig(2)) '}\)'])
zlabel(['\(f(x_{' int2str(whCoordBig(1)) '}, x_{' int2str(whCoordBig(2)) ...
   '}, ' num2str(nom_Val) ', \ldots)\)'])

figure
rotate3d on
f_app_Vis = reshape(f_app_Vis,[n_grd_val n_grd_val]);
surf(x_grd,x_grd,f_app_Vis); shading interp
zlim(zl)
xlabel(['\(x_{' int2str(whCoordBig(1)) '}\)'])
ylabel(['\(x_{' int2str(whCoordBig(2)) '}\)'])
zlabel(['\(f_{\mbox{app}}(x_{' int2str(whCoordBig(1)) '}, x_{' int2str(whCoordBig(2)) ...
   '}, ' num2str(nom_Val) ', \ldots)\)'])

figure
rotate3d on
err_Vis = f_true_Vis - f_app_Vis;
surf(x_grd,x_grd,err_Vis); shading interp
xlabel(['\(x_{' int2str(whCoordBig(1)) '}\)'])
ylabel(['\(x_{' int2str(whCoordBig(2)) '}\)'])
zlabel(['\(\mbox{err}(x_{' int2str(whCoordBig(1)) '}, x_{' int2str(whCoordBig(2)) ...
   '}, ' num2str(nom_Val) ', \ldots)\)'])

figure(1)

