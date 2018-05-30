%jac_h  measurement model of oinhole camera model with distortion coefficientts.
%
% Syntax:
%    J=jac_h(m,inde)
%
% In:
%   m  - Nx1 state 
%   inde - indexes of variables in state.
%  
% Out:
%   J  - jacobian matrix of pinhole camera projection
%
% Description:
%   Calculate the jacobian of the process of projecting the points into
%   the camera plane and apply distortion correction.
%
% Copyright (C) 2018 Santiago Cortés
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.


function J=jac_h(m,inde)
%Calculate jacobian matrix of measurement function h at point m.

%-dhp is the derivative of the position in the image plane (non
%normalized).
%-dh is de derivative of the pixel position before applying radial
%distortion.
%-dhr ir the derivative of the pixel position after radial distortion.



%Align vector
m = m(:);

%Number of points in state
inl_count=(length(m)-inde.w(end))/3;

%read, normalize points
X=m(inde.z(1):end);
X=reshape(X,[3,length(X)/3]);
Z=[X;ones(1,size(X,2))];

%extract rotation, position.
q=m(inde.q);
p=m(inde.p);

%Form intrinsic and extrinsic camera matrices.
f=m(inde.c);
K=[f(1) 0 f(3);0 f(2) f(4); 0 0 1];
R=quat2rmat(q');
E=[R' -R'*p];
Y=K*E*Z;

%derivative of image measurement with respect to camera parameters
dhp_df1=[1 0 0;0 0 0;0 0 0]*E*Z;
dhp_df2=[0 0 0;0 1 0;0 0 0]*E*Z;
dhp_dc1=[0 0 1;0 0 0;0 0 0]*E*Z;
dhp_dc2=[0 0 0;0 0 1;0 0 0]*E*Z;

%derivative with respect to 3d position of points and camera
dhp_dx=K*R';
dhp_dt=-K*R';
I=eye(3);

%pixel points before radial distortion.
im_point=(Y(1:2,:)./repmat(Y(3,:),2,1))';

%square radius
rad_sq=sum(((im_point-f(3:4)')./f(1:2)').^2,2);

%Derivative with respect to distortion coefficients.
dh_dk1=im_point.*rad_sq;
dh_dk2=im_point.*rad_sq.^2;
if inl_count==0
    J=[];
    return;
end
for i=1:inl_count
    
    %current point.
    x=X(:,i);
    
    %cross matrices.
    cross_x=[0 -x(3) x(2);x(3) 0 -x(1);-x(2) x(1) 0];
    cross_p=[0 -p(3) p(2);p(3) 0 -p(1);-p(2) p(1) 0];
    
    %Current pixel measurement.
    y=Y(:,i);
    y_pix=y(1:2)/y(3);
    
    %Current radius.
    r=rad_sq(i);
    
    v=-q(2:4);
    w=q(1);
    
    %Derivative of image measurement with respect to quaternion
    dhp_dq=2*K*[w*x+cross(v,x), (v'*x*I+v*x'-x*v'-w*cross_x)]-2*K*[w*p+cross(v,p), (v'*p*I+v*p'-p*v'-w*cross_p)];
    
    %Derivative of pixel measurement with respect to camera matrix
    %parameters
    dh_df1=dhp_df1(:,i)/y(3)-y/(y(3)^2)*dhp_df1(3,i);
    dh_df2=dhp_df2(:,i)/y(3)-y/(y(3)^2)*dhp_df2(3,i);
    dh_dc1=dhp_dc1(:,i)/y(3)-y/(y(3)^2)*dhp_dc1(3,i);
    dh_dc2=dhp_dc2(:,i)/y(3)-y/(y(3)^2)*dhp_dc2(3,i);
    dh_dc=[dh_df1, dh_df2, dh_dc1, dh_dc2];
    
    %Derivative of pixel measurement with respect to camera position and speed
    dh_dt=dhp_dt/y(3)-y/(y(3)^2)*dhp_dt(3,:);
    dh_dv=zeros(3);
    
    %apply cahin rule to derivatives to include normalization
    dh_dq=dhp_dq/y(3)-y/(y(3)^2)*dhp_dq(3,:);
    dh_dq(:,2:4)=-dh_dq(:,2:4);
    dh_dw=zeros(3);
    dh_dx=dhp_dx/y(3)-y/(y(3)^2)*dhp_dx(3,:);
    
    %Derivatives of radius measurement.
    dr_dc(1)=2*(y_pix(1)-f(3))*(dh_dc(1,1)*f(1)-y_pix(1))/f(1)^3    +2*(y_pix(2)-f(4))*(dh_dc(2,1))/f(2)^2;
    dr_dc(2)=2*(y_pix(1)-f(3))*(dh_dc(1,2))/f(1)^2                  +2*(y_pix(2)-f(4))*(dh_dc(2,2)*f(2)-y_pix(2))/f(2)^3;
    dr_dc(3)=2*(y_pix(1)-f(3))*(dh_dc(1,3)-1)/f(1)^2                +2*(y_pix(2)-f(4))*(dh_dc(2,3))/f(2)^2;
    dr_dc(4)=2*(y_pix(1)-f(3))*(dh_dc(1,4))/f(1)^2                  +2*(y_pix(2)-f(4))*(dh_dc(2,4)-1)/f(2)^2;
    dr_dt=2*(y_pix(1)-f(3))*dh_dt(1,:)/f(1)^2+2*(y_pix(2)-f(4))*dh_dt(2,:)/f(2)^2;
    dr_dv=2*(y_pix(1)-f(3))*dh_dv(1,:)/f(1)^2+2*(y_pix(2)-f(4))*dh_dv(2,:)/f(2)^2;
    dr_dq=2*(y_pix(1)-f(3))*dh_dq(1,:)/f(1)^2+2*(y_pix(2)-f(4))*dh_dq(2,:)/f(2)^2;
    dr_dw=2*(y_pix(1)-f(3))*dh_dw(1,:)/f(1)^2+2*(y_pix(2)-f(4))*dh_dw(2,:)/f(2)^2;
    dr_dx=2*(y_pix(1)-f(3))*dh_dx(1,:)/f(1)^2+2*(y_pix(2)-f(4))*dh_dx(2,:)/f(2)^2;
    
    
    %Apply cahin rule to derivatives to include radial distortion.
    for dim=1:2
        dh_dcr(dim,:)=dh_dc(dim,:)*(1+f(5)*r+f(6)*r^2)+y_pix(dim)*dr_dc*(f(5)+f(6)*2*r);
        dh_dtr(dim,:)=dh_dt(dim,:)*(1+f(5)*r+f(6)*r^2)+y_pix(dim)*dr_dt*(f(5)+f(6)*2*r);
        dh_dvr(dim,:)=dh_dv(dim,:)*(1+f(5)*r+f(6)*r^2)+y_pix(dim)*dr_dv*(f(5)+f(6)*2*r);
        dh_dqr(dim,:)=dh_dq(dim,:)*(1+f(5)*r+f(6)*r^2)+y_pix(dim)*dr_dq*(f(5)+f(6)*2*r);
        dh_dwr(dim,:)=dh_dw(dim,:)*(1+f(5)*r+f(6)*r^2)+y_pix(dim)*dr_dw*(f(5)+f(6)*2*r);
        dh_dxr(dim,:)=dh_dx(dim,:)*(1+f(5)*r+f(6)*r^2)+y_pix(dim)*dr_dx*(f(5)+f(6)*2*r);
    end
    
    %Put results in matrix.
    index_in_matrix=[i,i+inl_count];
    J(index_in_matrix,inde.c(1:4))=dh_dcr(1:2,:);
    J(index_in_matrix,inde.c(5:6))=[dh_dk1(i,:)', dh_dk2(i,:)'];
    J(index_in_matrix,inde.p)=dh_dtr(1:2,:);
    J(index_in_matrix,inde.v)=dh_dvr(1:2,:);
    J(index_in_matrix,inde.q)=dh_dqr(1:2,:);
    J(index_in_matrix,inde.w)=dh_dwr(1:2,:);
    J(index_in_matrix,inde.z(1):inde.z(1)+inl_count*3-1)=0;
    J(index_in_matrix,inde.z(i*3-2):inde.z(i*3-2)+2)=dh_dxr(1:2,:);
    
end
return

%% Check that derivatives match
clear all
rng(0,'twister')
%clc

% Indices
inde.c=1:6;
inde.p=inde.c(end)+1:inde.c(end)+3;
inde.v=inde.p(end)+1:inde.p(end)+3;
inde.pv=[inde.p inde.v];
inde.q=inde.v(end)+1:inde.v(end)+4;
inde.w=inde.q(end)+1:inde.q(end)+3;
N_feat=8;
inde.z=inde.w(end)+1:inde.w(end)+N_feat*3;

% Give the parameters some values

init(3,:)=randn([N_feat,1])-10;
init(1,:)=randn([N_feat,1]);
init(2,:)=randn([N_feat,1]);
%init=init';
x = randn(4+6+4+N_feat*3+3,1);


%point camera to one of the features
dx = init(1,1);
dy = init(2,1);
dz = init(3,1);

% Define vectors
a=[0 0 1];
b=-[dx dy dz];
b=b/norm(b);

% Find axis of rotation
r_axis=cross(a,b);
r_axis=r_axis/norm(r_axis);

% Angle between vectors
%r_angle = atan2(norm(cross(a,b)),dot(a,b));
r_angle = acos(dot(a,b));

% Compute correspoding quaternion
q(1)=cos(r_angle/2);
q(2)=r_axis(1)*sin(r_angle/2);
q(3)=r_axis(2)*sin(r_angle/2);
q(4)=r_axis(3)*sin(r_angle/2);
q=q/norm(q);




x(inde.q) = q;

x(inde.c) = 500;
x(inde.z)=init(:);
x(inde.p)=0;
x(inde.c(5))=1e-3;
x(inde.c(6))=1e-6;
% Check derivatives
[D0,D1] = der_check(@(x) h(x,inde), ...
    @(x) jac_h(x,inde),1,x);

der_check(@(x) h(x,inde), ...
    @(x) jac_h(x,inde),1,x);

sum(sum(abs(D0-D1)))

figure(1); clf
imagesc(log10(abs(D0-D1)))
axis equal tight
colorbar
of=1;
test=inde.c