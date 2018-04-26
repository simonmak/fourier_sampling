function [basisVal] = eval_f_four_add(xnew,basisVal,basisFun,waveNum,s_max)
% Adds polynomial evaluations for a new design point xnew
nBasis = size(waveNum,1);
% [nX,d] = size(x);
d = length(xnew);
[nX,~,~] = size(basisVal); %add +1 to design size
nX = nX + 1;

newSlice = zeros(1,d,s_max+1); %initialize new slice to add to basisVal
newSlice(1,:,1) = 1; %first value is always the constant one
p_vec = ones(1,s_max+1);
for s = 1:s_max
    [newSlice(1,:,s+1),p_vec(s+1)] = basisFun(s,xnew); %evaluate basis functions at xnew
end
% pVal(nBasis,1) = 0;
% for j = 1:nBasis
%     pVal(j) = prod(p_vec(waveNum(j,:)+1));
% end

basisVal = cat(1,basisVal,newSlice); %add new slice to polynomial evaluations