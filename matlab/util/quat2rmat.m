function R = quat2rmat( q )
% Converts quaternions into rotation matrices
% If a vector of quaternions is given, it needs to be stacked vertically
% Uses a different implementation for one or multiple quaternions for
% computational efficiency

if any(size(q) == 1)
    q0 = q(1); q1 = q(2); q2 = q(3); q3 = q(4);

    R = [q0^2+q1^2-q2^2-q3^2, 2*q1*q2 - 2*q0*q3 ,2*q1*q3 + 2*q0*q2 ; ...
        2*q1*q2 + 2*q0*q3, q0^2 - q1^2 + q2^2 - q3^2, 2*q2*q3 - 2*q0*q1 ; ...
       2*q1*q3 - 2*q0*q2, 2*q2*q3 + 2*q0*q1, q0^2 - q1^2 - q2^2 + q3^2];
else
    nQuat = size(q,1);
    qvec = reshape(q(:,2:4)',[3,1,nQuat]);
    R = multiprod(qvec, reshape(q(:,2:4)',[1,3,nQuat])) + ...
        multiprod(reshape(q(:,1).^2,1,1,nQuat),eye(3)) + ...
        2*multiprod(reshape(q(:,1),1,1,nQuat),mcross(q(:,2:4))) + ...
        multiprod(mcross(q(:,2:4)),mcross(q(:,2:4)));
end

end