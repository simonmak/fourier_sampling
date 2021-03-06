function [] = sim_direct_plot_results(d,s_flg)

    %%Plot results
%     clearvars

    load(['sim_eval_results_d' int2str(d) '_sflg' int2str(s_flg) '.mat']) %load results for plotting
    InitializeDisplay %add some variables for nice plotting

    figure
    log10epsVec = log10(eps_vec);
    h = scatter(rat_vec,n_vec,800,log10epsVec,'.'); %plot ratio of actual error to tolerance, with color corresonding to tolerance
    hold on
    h = [h; scatter(rat_four_eps_vec,n_vec,300,log10epsVec,'x','linewidth',4)]; %ratio of L^1 error of four_coef to tolerance
    h = [h; scatter(rat_normf_eps_vec,n_vec,300,log10epsVec,'+','linewidth',4)]; %ratio of L^1 error of gamma to tolerance
    title( ['Random Fourier \((d = ' num2str(d) ')\)'] )
    set(gca,'XScale','log', 'YScale','log')
    xlims = [min([[rat_vec; rat_four_eps_vec; rat_normf_eps_vec]*0.8; 0.1]) ...
     max([[rat_vec; rat_four_eps_vec; rat_normf_eps_vec]*1.2; 1])];
    xlim(xlims)
    disp(xlims)
    xlim([0.0067 1.1998])
    ylims = 10.^[floor(log10(min(n_vec))) ceil(log10(max(n_vec)))];
    ylim(ylims)
    disp(ylims)
    ylim([10 100000])
    xlabel({'\(||f-\hat{f}||_{\infty}/\varepsilon\)'})
    ylabel({'Sample size \(n\)'})
    hcb = colorbar; %showing tolerance values
    title(hcb,'\(\varepsilon\)','interpreter','latex')
    tickVals = floor(min_log10_eps:max_log10_eps);
    tickLabels = 10.^tickVals;
    set(hcb,'Ticks',tickVals,'TickLabels',tickLabels, ...
       'Limits',[tickVals(1) tickVals(end)])
    [~,leg_icons] = legend(h,{'\(||f-f_{\mbox{app}}||_{\infty}/\varepsilon\) \quad', ...
       '\(||\hat{f} - \hat{f}_{\mbox{app}}||_{1}/\varepsilon\)\quad', ...
       '\(||\hat{f}||_{\infty,}\)\boldmath\({}_\gamma\)\unboldmath\(||\bigl(\gamma(\)\boldmath\(j\)\unboldmath\({}_i) \bigr)_{i > n}||_{1}/\varepsilon\)'}, ...
       'box','off','location','south','orientation','horizontal');
    set(gcf,'Position',[200,200,1000,500]) %make figure big enough and the right aspect ratio
    leg_icons(4).Children.MarkerSize = 30; %make legend icons large enough
    leg_icons(5).Children.MarkerSize = 20;
    leg_icons(6).Children.MarkerSize = 20;
    leg_icons(5).Children.LineWidth = 4;
    leg_icons(6).Children.LineWidth = 4;
    %legend boxoff

    %% Visualize (only a 2-d projection)
    [~,whCoordBig] = sort(w_vec);
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
    [f_true_Vis,basisVal] = eval_f_four(x_Vis,basisFun,gam_mtx,s_max,four_coef);
    [f_app_Vis] = ...
           eval_f_four([],basisVal,gam_mtx(gam_idx(1:nn),:),s_max,four_coef(gam_idx(1:nn)));

    figure
    rotate3d on
    f_true_Vis = reshape(f_true_Vis,[n_grd_val n_grd_val]);
    surf(x_grd,x_grd,f_true_Vis); shading interp
    xlabel(['\(x_{' int2str(whCoordBig(1)) '}\)'])
    ylabel(['\(x_{' int2str(whCoordBig(2)) '}\)'])
    zlabel(['\(f(x_{' int2str(whCoordBig(1)) '}, x_{' int2str(whCoordBig(2)) ...
       '}, ' num2str(nom_Val) ', \ldots)\)'])

    figure
    rotate3d on
    f_app_Vis = reshape(f_app_Vis,[n_grd_val n_grd_val]);
    surf(x_grd,x_grd,f_app_Vis); shading interp
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

end