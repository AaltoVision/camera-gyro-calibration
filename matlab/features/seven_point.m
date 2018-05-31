%seven_point  Seven point algorithm for fundamental matrix
%
% Syntax:
%    S=seven_point(P)
%
% In:
%   P  - Points
%   
%  
% Out:
%   S  - fundamental matrix
%
% Description:
%   Seven point algorithm for estimation of a fundamental matrix qiven two
%   sets of corresponding 7 points. Following Hartley Zisserman book.
%
%   Richard Hartley and Andrew Zisserman. 2003. Multiple View Geometry 
%       in Computer Vision (2 ed.). Cambridge University Press, New York,
%       NY, USA
%
% Copyright (C) 2018 Santiago Cortés
%
% This software is distributed under the GNU General Public 
% Licence (version 2 or later); please refer to the file 
% Licence.txt, included with the software, for details.

function S=seven_point(P)


%form matrix
for i=1:7
    v=([P(i,3:4) 1]'*[[P(i,1:2) 1]])';
    A(i,:)=v(:) ;
end
%find left null space of matrix
N=null(A);

%generate F matrices
F1=reshape(N(:,1),[3 3])';
F2=reshape(N(:,2),[3 3])';
%solve for te correct fundamental mmatrix
x=sym('x');
y=det(x*F1+(1-x)*F2);
c=sym2poly(y);
R=roots(c);
S=[];
ind=1;
for i=1:3
    if abs(imag(R(i)))<1e-16
    S(:,:,ind)=(real(R(i))*F1+(1-real(R(i)))*F2);
    ind=ind+1;
    else
       1; 
    end
end
1;
