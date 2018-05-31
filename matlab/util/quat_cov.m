%quat_cov process noise for quaternion given a rotational rate. 
%
% Syntax:
%   S= quat_cov(q,w_var,del_t)
%
% In:
%   q  - unit quaternion representing initial orientation.
%   w_var - variance of the rotational rate (confidence on measurement)
%   del_t - time window (constant rotation rate assumed)
%
%  
% Out:
%   S  - covariance matrix
%
% Description:
%   Propagate uncertainty from rotational rate into rotation quaternion.
%
% Copyright (C) 2018 Santiago Cortés
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.
function S= quat_cov(q,w_var,del_t)
%process noise for quaternion given a rotational rate.
M=[-q(2:4);q(1) -q(4) q(3);q(4) q(1) -q(2);-q(3) q(2) q(1)];
S=del_t^2*1/4*M*diag(w_var)*M';


