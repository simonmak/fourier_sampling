function [basisVal,XX] = eval_X_add(xnew,basisVal,XX,basisFun,waveNum,s_max,nPts,nBasisDes)
% Adds polynomial evaluations for a new design point xnew
d = size(xnew,2);

%% Update basisVal
newSlice = zeros(1,d,s_max+1); %initialize new slice to add to basisVal
newSlice(1,:,1) = 1; %first value is always the constant one
%p_vec = ones(1,s_max+1);
for s = 1:s_max
    newSlice(1,:,s+1) = basisFun(s,xnew); %evaluate basis functions at xnew
end
basisVal(nPts,:,:) = newSlice; %add new slice to polynomial evaluations

%% Update design matrix

%Update new point for all bases
for j = 1:nBasisDes
   addPart = 1;
   for ell = 1:d
      addPart = addPart .* basisVal(nPts,ell,waveNum(j,ell)+1);
   end
   XX(nPts,j) = addPart;
end

%Update new basis for old points
addPart = 1;
for ell = 1:d
   addPart = addPart .* basisVal(1:nPts-1,ell,waveNum(nBasisDes,ell)+1);
end
XX(1:nPts-1,nBasisDes) = addPart;

