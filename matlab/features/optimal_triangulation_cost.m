%optimal_triangulation_cost  Optimal cost for triangulation ZIzzerman
%algorithm 12.1
%
% Syntax:
%   S=optimal_triangulation_cost(x,xp,F)
%
% In:
%   x  - pixel coordinates in camera 1 
%   xp  - pixel coordinates in camera 2
%   F - potential Fundamental matrix
%  
% Out:
%   S  - cost to be minimized
%
% Description:
%   Part of algorithm 12.2 in the book:
%   Richard Hartley and Andrew Zisserman. 2003. Multiple View Geometry 
%       in Computer Vision (2 ed.). Cambridge University Press, New York,
%       NY, USA
%
% Copyright (C) 2018 Santiago Cortés
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function S=optimal_triangulation_cost(x,xp,F)
%cost of match given a fundamental matrix, Zisserman book algorithm 12.1
T=eye(3);
Tp=eye(3);
T(1:2,3)=-x;
Tp(1:2,3)=-xp;
F=(Tp')^-1*F*T^-1;
e=null(F);
ep=null(F');
if isempty(e) || isempty(ep)
    S=inf;
    return 
end
e=e/norm(e(1:2));
ep=ep/norm(ep(1:2));
R=[e(1) e(2) 0;-e(2) e(1) 0;0 0 1];
Rp=[ep(1) ep(2) 0;-ep(2) ep(1) 0;0 0 1];
F=Rp*F*R';
f=e(3);
v=ep(3);
a=F(2,2);
b=F(2,3);
c=F(3,2);
d=F(3,3);
g(1)=   a*b*c^2*f^4     -   f^4*d*c*a^2;
g(2)=   c^4*v^4        +   2*a^2*c*d*v^2  +   4*a*b*c^2*v^2  -   a*b*d^2*f^4 ...
    +   b^2*c*d*f^4     -   2*a^2*c*d*f^2   +   2*a*b*c^2*f^2   +   4*a^3*b;
g(3)=   4*c^3*d*v^4     +   4*a^2*c*d*v^2   +   4*a*b*c^2*v^2   -   a*b*d^2*f^4 ...
    +   b^2*c*d*f^4     -   2*a^2*c*d*f^2   +   2*a*b*c^2*f^2   +   4*a^3*b;
g(4)=   6*c^2*d^2*v^4   +   2*a^2*d^2*v^2   +   8*a*b*c*d*v^2   +   2*b^2*c^2*v^2  ...
    -   2*a^2*d^2*f^2   +   2*b^2*c^2*f^2   +   6*a^2*b^2;
g(5)=   4*c*d^3*v^4     +   4*a*b*d^2*v^2   +   4*b^2*c*d*v^2   -   2*a*b*d^2*f^2   ...
    +   2*b^2*c*d*f^2   -   a^2*c*d         +   a*b*c^2         +   4*a*b^3;
g(6)=   d^4*v^4         +   2*b^2*d^2*v^2   -   a^2*d^2+b^2*c^2 +   b^4;
g(7)=   b^2*c*d         -   a*b*d^2;
tv=roots(g);

for i=1:6
    t=real(tv(i));
    
    s(i)=       (  (t^2)   /   (1+(f^2)*(t^2)) )    +   ...
        (   ((c*t+d)^2) /   (   ((a*t+b)^2) +   (v^2)*((c*t+d)^2)   )    );
    
end
%s(7)=(1/(f^2))+(c^2)/((a^2)+(v^2)*(c^2));
[~,ind]=min(s);
tmin=tv(ind);
l=[tmin*f,1,-tmin]';
lp=[-v*(c*t+d),a*t+b,c*t+d]';

xt=[-l(1)*l(3),-l(2)*l(3),l(1)^2+l(2)^2]';
xtp=[-lp(1)*lp(3),-lp(2)*lp(3),lp(1)^2+lp(2)^2]';

xt=(T)^-1*R'*xt;
xtp=(Tp)^-1*Rp'*xtp;

xt=xt(1:2)/xt(3);
xtp=xtp(1:2)/xtp(3);

S=norm(x-xt')+norm(xtp-xp');
1;

