function [basisVal,XX] = eval_X_add(xnew,basisVal,XX,basisFun,waveNum,s_max)
% Adds polynomial evaluations for a new design point xnew
d = size(xnew,2);
[nPts,nBasisDes] = size(XX);

%% Update basisVal
newSlice = zeros(1,d,s_max+1); %initialize new slice to add to basisVal
newSlice(1,:,1) = 1; %first value is always the constant one
%p_vec = ones(1,s_max+1);
for s = 1:s_max
    newSlice(1,:,s+1) = basisFun(s,xnew); %evaluate basis functions at xnew
end
basisVal = cat(1,basisVal,newSlice); %add new slice to polynomial evaluations

%% Update design matrix
XX = [XX, zeros(nPts,1); zeros(1, nBasisDes+1)];

%Update new point for all bases
for j = 1:nBasisDes
   addPart = 1;
   for ell = 1:d
      addPart = addPart .* basisVal(nPts+1,ell,waveNum(j,ell)+1);
   end
   XX(nPts+1,j) = addPart;
end

%Update new basis for old points
addPart = 1;
for ell = 1:d
   addPart = addPart .* basisVal(:,ell,waveNum(nBasisDes+1,ell)+1);
end
XX(:,nBasisDes+1) = addPart;

