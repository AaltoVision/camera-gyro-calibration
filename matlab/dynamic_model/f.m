%f  dynamic model of inertial motion system asuming w as a control signal.
%
% Syntax:
%   S = f(x,delT,inde)
%
% In:
%   x  - Nx1 state 
%   delT  - 1x1 time to integrate.
%   inde - indexes of variables in state.
%  
% Out:
%   S  - Updated state mean
%
% Description:
%   Propagate the state using the wiener velocity model.
%
% Copyright (C) 2018 Santiago Cortés
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function S = f(x,delT,inde)
%c=fx,fy,cx,cy.
%p is position in R3.
%q is rotation as a quaternion.
%w is the gyroscope data, in R3.
%z is the feature vector positions (N by 3).
%x= c p q z;

qo=x(inde.q);
w=x(inde.w);
qo=qo/norm(qo);
%qn=circProd(qo,w,delT);
crow=[0 w(3) -w(2);-w(3) 0 w(1);w(2) -w(1) 0];
ome=[0 -w;w' crow ];
qn=(expm(delT*ome/2)*qo')';
S=x;
S(inde.q)=qn;
S(inde.p)=S(inde.p)+S(inde.v)*delT;
