%circProd  apply rotation of rotational speed w over time delT to quaternion q.
%
% Syntax:
%   S=circProd(q,w,delT)
%
% In:
%   q  - 4x1 initial orientation quaternion.
%   w  - 3x1 rotational rate measurement.
%   delT  - 1x1 time to integrate.
%  
% Out:
%   s  - Updated state mean
%
% Description:
%   Assume constant rotational rate over period of delT time and calculate
%   resulting orientation quaternion
%


% Copyright (C) 2018 Santiago Cortés
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function S=circProd(q,w,delT)
%apply rotation of rotational speed w over time delT to quaternion q.

%avoid operation for null rotation
if norm(w)==0
    S=q;
else
    %zero order forward integrator for rotation rate in quaternions.
    ys=w/norm(w)*sin(norm(w)*delT/2);
    yc=cos(norm(w)*delT/2);
    S(1)=q(1)*yc-dot(q(2:4),ys);
    S(2:4)=q(1)*ys+yc*q(2:4)+cross(ys,q(2:4));
end