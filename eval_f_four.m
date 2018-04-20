function [fVal,basisVal,pVal] = eval_f_four(x,basisVal,waveNum,s_max,fourCoeff)
% Evaluate a function based on Fourier coefficient values
nBasis = size(waveNum,1);
if isa(basisVal,'function_handle') %Don't yet have basis tablulated at x
   [nX,d] = size(x);
   basisFun = basisVal; %the values of the basis is really the basis function
   basisVal = zeros(nX,d,s_max+1); %initialize the values
   basisVal(:,:,1) = 1; %first value is always the constant one
   p_vec = ones(1,s_max+1);
   for s = 1:s_max
     [basisVal(:,:,s+1),p_vec(s+1)] = basisFun(s,x); %evaluate basis functions at x
   end
   pVal(nBasis,1) = 0;
   for j = 1:nBasis
      pVal(j) = prod(p_vec(waveNum(j,:)+1));
   end
else
   [nX,d,~] = size(basisVal);
   pVal = [];
end

fVal(nX,1) = 0;
for j = 1:nBasis
   addPart = ones(nX,1);
   for ell = 1:d
      addPart = addPart .* basisVal(:,ell,waveNum(j,ell)+1);
   end
   fVal = fVal + fourCoeff(j)*addPart;
end