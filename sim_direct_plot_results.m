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
tickVals = floor(min(log10epsVec)):ceil(max(log10epsVec));
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
lp_Vis(n_Vis,d,s_max+1) = 0;
lp_Vis(:,:,1) = 1;
for s = 1:s_max
  temp = legendre(s,x_Vis); %generate associated Legendre functions
  lp_Vis(:,:,s+1) = squeeze(temp(1,:,:)); %keep Legendre polynomials
end

%Evaluate f_true
f_true_Vis(n_Vis,1) = 0;
for j = 1:nBasis
   addPart = ones(n_Vis,1);
   for ell = 1:d
      addPart = addPart .* lp_Vis(:,ell,gam_mtx(j,ell)+1);
   end
   f_true_Vis = f_true_Vis + four_coef(j)*addPart;
end

%Evaluate f_app
f_app_Vis = zeros(n_Vis,1);
for j = 1:nn
 jj = gam_idx(j); %which basis is next smallest
 addPart = ones(n_Vis,1);
 for ell = 1:d
    addPart = addPart .* lp_Vis(:,ell,gam_mtx(jj,ell)+1);
 end
 f_app_Vis = f_app_Vis + four_coef(jj)*addPart;
end


figure
rotate3d on
f_true_Vis = reshape(f_true_Vis,[n_grd_val n_grd_val]);
surf(x_grd,x_grd,f_true_Vis); shading interp

figure
rotate3d on
f_app_Vis = reshape(f_app_Vis,[n_grd_val n_grd_val]);
surf(x_grd,x_grd,f_app_Vis); shading interp

figure(1)

