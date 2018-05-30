%h  measurement model of oinhole camera model with distortion coefficientts.
%
% Syntax:
%    S = h(x,inde)
%
% In:
%   x  - Nx1 state 
%   inde - indexes of variables in state.
%  
% Out:
%   S  - pixel measurements
%
% Description:
%   Project the points into the camera plane and apply distortion
%   correction
%
% Copyright (C) 2018 Santiago Cortés
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.


function S = h(x,inde)
%Measure the pixel position of the points.

%c=fx,fy,cx,cy,k1,k2.
%p is position in R3.
%q is rotation as a quaternion.
%w is the gyroscope data, in R3.
%z is the feature vector positions (N by 3).
%x= c p q z;

%number of points in state.
inl_count=(size(x,1)-inde.w(end));

%if several states are presented, evaluate and return several measurments.
if inl_count==0
    S=[];
    return
end
if size(x,2)>1
    S = nan(2*(inl_count/3),size(x,2));
    for i=1:size(x,2)
        S(:,i) = h(x(:,i),inde);
    end
    return
end

%Align vector.
x=x(:);

%Extract positions and form homogenous matrix
z=x(inde.z(1):end);
Z=reshape(z,[3,length(z)/3]);
Z=[Z;ones(1,size(Z,2))];
%Extract and form intrinsic camera matrix
f=x(inde.c);
K=[f(1) 0 f(3);0 f(2) f(4); 0 0 1];

%Extract and normalize quaternion, then form matrix.
q=x(inde.q)/norm(x(inde.q));
R=quat2rmat(q);

%Extract camera position and form extrinsic matrix.
p=x(inde.p);
E=[R' -R'*p];

%Project onto camera plane
Y=K*E*Z;

%Normalize to pixel measurement
S=(Y(1:2,:)./repmat(Y(3,:),2,1))';

%calculate radius in normalized units.
r=sum(((S-f(3:4)')./f(1:2)').^2,2);

%apply radial distortion.
S(:,1)=S(:,1).*(1+f(5)*r+f(6)*r.^2);
S(:,2)=S(:,2).*(1+f(5)*r+f(6)*r.^2);

%Form measurement vector.
S=S(:);
