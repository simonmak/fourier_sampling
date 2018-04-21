function [y,p] = chebyshevBasis(s,x)
y = cos(s*acos(x));
p = 1;
return
