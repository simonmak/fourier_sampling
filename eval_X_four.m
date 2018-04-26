function [XX,basisVal,pVal] = eval_X_four(basisVal,waveNum)
% Evaluate model matrix X based on pre-computed basisVal
nBasis = size(waveNum,1);
[nX,d,~] = size(basisVal);
pVal = [];

XX = zeros(nX,nBasis);
for j = 1:nBasis
   addPart = ones(nX,1);
   for ell = 1:d
      addPart = addPart .* basisVal(:,ell,waveNum(j,ell)+1);
   end
   XX(:,j) = addPart;
end