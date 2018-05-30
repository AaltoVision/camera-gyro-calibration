function S= quat_cov(q,w_var,del_t)
%process noise for quaternion given a rotational rate.
M=[-q(2:4);q(1) -q(4) q(3);q(4) q(1) -q(2);-q(3) q(2) q(1)];
S=del_t^2*1/4*M*diag(w_var)*M';



%%
% for i=1:4
%    for j=1:4
%      S(i,j)=-1/12*q_cov*del_t^3*q(i)*q(j);  
%    end
% end
% S(1,1)=1/12*q_cov*del_t^3*(sum(q([4 3 2]).^2));
% S(2,2)=1/12*q_cov*del_t^3*(sum(q([4 3 1]).^2));
% S(3,3)=1/12*q_cov*del_t^3*(sum(q([4 2 1]).^2));
% S(4,4)=1/12*q_cov*del_t^3*(sum(q([3 2 1]).^2));
% 1;

