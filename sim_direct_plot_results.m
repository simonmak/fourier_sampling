%%Plot results
clearvars

load sim_direct_results.mat %load results for plotting
InitializeDisplay %add some variables for nice plotting

figure
log10epsVec = log10(eps_vec);
scatter(rat_vec,n_vec,800,log10epsVec,'.')
set(gca,'YScale','log')
xlim([0 1])
ylim(10.^[floor(log10(min(n_vec))) ceil(log10(max(n_vec)))])
xlabel({'\(||f-\hat{f}||_{\infty}/\epsilon\)'})
ylabel({'Sample size \(n\)'})
hcb = colorbar;
title(hcb,'\(\epsilon\)','Interpreter','latex')
tickVals = floor(min(log10epsVec)):ceil(max(log10epsVec));
tickLabels = 10.^tickVals;
set(hcb,'Ticks',tickVals,'TickLabels',tickLabels, ...
   'Limits',[tickVals(1) tickVals(end)])

% %Visualize (only in 2-d)
% n_grd = 100;
% x_grd = -1:(2/n_grd):1;
% f_true = zeros(length(x_grd),length(x_grd));
% f_app = zeros(length(x_grd),length(x_grd));
% 
% for (i = 1:size(gam_mtx,1))
%     cur_idx = gam_mtx(i,:);
%     run_sum = gam_val(i)*(legendreP(cur_idx(1),x_grd))'*legendreP(cur_idx(2),x_grd);
%     %run_sum = run_sum / sqrt(2/(2*cur_idx(1)+1)) / sqrt(2/(2*cur_idx(2)+1)); %normalize
%     f_true = f_true + run_sum;
%     if (ismember(i,gam_idx(1:nn)))
%         f_app = f_app + run_sum;
%     end
% end
% 
% % figure
% % rotate3d on
% % [plotX,plotY] = meshgrid(x_grd);
% % mesh(plotX,plotY,f_true)
% % 
% % figure
% % rotate3d on
% % [plotX,plotY] = meshgrid(x_grd);
% % mesh(plotX,plotY,f_app)


