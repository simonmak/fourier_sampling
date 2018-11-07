function [basisVal,nX,d] = eval_Basis(x,basisFun,s_max)
% Evaluate a function based on Fourier coefficient values
[nX,d] = size(x);
basisVal = zeros(nX,d,s_max+1); %initialize the values
% basisVal = zeros(d,s_max+1,nX); %initialize the values
basisVal(:,:,1) = 1; %first value is always the constant one
% basisVal(:,1,:) = 1; %first value is always the constant one
for s = 1:s_max
  basisVal(:,:,s+1) = basisFun(s,x); %evaluate basis functions at x
%   basisVal(:,s+1,:) = basisFun(s,x)'; %evaluate basis functions at x
end
