function [fVal,basisVal] = eval_f_four(x,basisVal,waveNum,s_max,fourCoeff)
% Evaluate a function based on Fourier coefficient values
nBasis = size(waveNum,1);
if isa(basisVal,'function_handle') %Don't yet have basis tablulated at x
   basisFun = basisVal;
   [basisVal,nX,d] = eval_Basis(x,basisFun,s_max);
else
   [nX,d,~] = size(basisVal);
end

fVal(nX,1) = 0;
for j = 1:nBasis
   addPart = ones(nX,1);
   for ell = 1:d
      addPart = addPart .* basisVal(:,ell,waveNum(j,ell)+1);
   end
   fVal = fVal + fourCoeff(j)*addPart;
end