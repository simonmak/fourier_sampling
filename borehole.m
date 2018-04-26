function [y] = borehole(xx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% BOREHOLE FUNCTION
%
% Authors: Sonja Surjanovic, Simon Fraser University
%          Derek Bingham, Simon Fraser University
% Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
%
% Copyright 2013. Derek Bingham, Simon Fraser University.
%
% THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
% FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
% derivative works, such modified software should be clearly marked.
% Additionally, this program is free software; you can redistribute it 
% and/or modify it under the terms of the GNU General Public License as 
% published by the Free Software Foundation; version 2.0 of the License. 
% Accordingly, this program is distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
% of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
% General Public License for more details.
%
% For function details and reference information, see:
% http://www.sfu.ca/~ssurjano/
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% OUTPUT AND INPUT:
%
% y  = water flow rate
% xx = [rw, r, Tu, Hu, Tl, Hl, L, Kw]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rw = 0.05*xx(1)+0.1;
r  = 24950*xx(2)+25050;
Tu = 26265*xx(3)+89335;
Hu = 60*xx(4)+1050;
Tl = 26.45*xx(5)+89.55;
Hl = 60*xx(6)+760;
L  = 280*xx(7)+1400;
Kw = 1095*xx(8)+10950;

frac1 = 2 * pi * Tu * (Hu-Hl);

frac2a = 2*L*Tu / (log(r/rw)*rw^2*Kw);
frac2b = Tu / Tl;
frac2 = log(r/rw) * (1+frac2a+frac2b);

y = frac1 / frac2;

end
