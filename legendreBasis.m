function [y,p] = legendreBasis(s,x)
temp = legendre(s,x); %generate associated Legendre functions
y = squeeze(temp(1,:,:)); %keep Legendre polynomials
p = 1;
return
