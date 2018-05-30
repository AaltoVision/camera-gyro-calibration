function S=seven_point(P)
%Seven point algorithm for estimation of a fundamental matrix qiven two
%sets of corresponding 7 points. Following Hartley Zisserman book.

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
