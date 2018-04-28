function [XX,basisVal] = eval_X_four(x,basisVal,waveNum,s_max)
% Evaluate model matrix X at points x
nBasis = size(waveNum,1);
if isa(basisVal,'function_handle') %Don't yet have basis tablulated at x
   basisFun = basisVal;
   [basisVal,nX,d] = eval_Basis(x,basisFun,s_max);
else
   [nX,d,~] = size(basisVal);
end

XX = zeros(nX,nBasis);
for j = 1:nBasis
   addPart = ones(nX,1);
   for ell = 1:d
      addPart = addPart .* basisVal(:,ell,waveNum(j,ell)+1);
   end
   XX(:,j) = addPart;
end

