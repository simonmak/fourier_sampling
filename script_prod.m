% Simple script testing product weight opt.

d = 5;     %dimension
Gam_vec = 1./factorial(0:d); %\Gamma (order wts)
s_max = 3; %maximum smoothness
s_vec = 1./(( (0:s_max) +1).^2); %s (smoothness wts)
eps = 0.1;

% Fourier coef. setting for true f'n f
Gam_vec_tr = 1./factorial(0:d); %\Gamma (order wts)
w_vec_tr = 1./((1:d).^2); %w (product wts)
s_max_tr = 3; %maximum smoothness
s_vec_tr = 1./(( (0:s_max_tr) +1).^2); %s (smoothness wts)

% Compute Fourier coef and gammas
rad_seq = ones(length(four_coef),1);
four_coef = comp_wts(Gam_vec_tr,w_vec_tr,s_vec_tr,s_max_tr);
four_coef = rad_seq .* four_coef;
C = 1.2;
n0 = 3;
w_ini = 0.25*ones(d,1); %init. for w optimization

% Estimate product weights
w_vec = -1*ones(d,1); %set dummy weights if no flag
[~,~,w_est] = samp_sz(four_coef,Gam_vec,w_vec,s_vec,s_max,[],eps,1.0,n0,false,false,w_ini);
disp(w_est) %estimated weights via nonlinear opt.
disp(w_vec_tr') %true weights

