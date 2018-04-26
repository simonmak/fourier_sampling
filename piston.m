function [C] = piston(xx)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% PISTON SIMULATION FUNCTION
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
% C = cycle time
% xx = [M, S, V0, k, P0, Ta, T0]
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M  = 15*xx(1)+45;
S  = 0.0075*xx(2)+0.0125;
V0 = 0.004*xx(3)+0.006;
k  = 2000*xx(4)+3000;
P0 = 10000*xx(5)+100000;
Ta = 3*xx(6)+293;
T0 = 10*xx(7)+350;

Aterm1 = P0 * S;
Aterm2 = 19.62 * M;
Aterm3 = -k*V0 / S;
A = Aterm1 + Aterm2 + Aterm3;

Vfact1 = S / (2*k);
Vfact2 = sqrt(A^2 + 4*k*(P0*V0/T0)*Ta);
V = Vfact1 * (Vfact2 - A);

fact1 = M;
fact2 = k + (S^2)*(P0*V0/T0)*(Ta/(V^2));

C = 2 * pi * sqrt(fact1/fact2);

end
