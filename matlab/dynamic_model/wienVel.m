%wienVel process noise for wiener velocity model given a rotational rate
%
% Syntax:
%    Q=wienVel(w,delT)
%
% In:
%   w  - rotational rate
%   delT - time window
%  
% Out:
%   Q  - process noise 
%
% Description:
%   Calculate process noise for orientation given a rotational rate
%
% Copyright (C) 2018 Santiago Cortés
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function Q=wienVel(w,delT)
%Wiener velocity model process noise for given roation rate.
if length(w==1)
    w=[w w w];
end
Q=[diag(w*(delT^3)/3) diag(w*(delT^2)/2);diag(w*(delT^2)/2) diag(w*delT)];