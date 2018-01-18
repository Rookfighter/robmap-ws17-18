% Compute the error of a pose-pose constraint
% x1 3x1 vector (x,y,theta) of the first robot pose
% x2 3x1 vector (x,y,theta) of the second robot pose
% z 3x1 vector (x,y,theta) of the measurement
%
% You may use the functions v2t() and t2v() to compute
% a Homogeneous matrix out of a (x, y, theta) vector
% for computing the error.
%
% Output
% e 3x1 error of the constraint
% A 3x3 Jacobian wrt x1
% B 3x3 Jacobian wrt x2
function [e, A, B] = linearize_pose_pose_constraint(x1, x2, z)

    % transform poses and meaurement to homog corrds
    xt1 = v2t(x1);
    xt2 = v2t(x2);
    zt  = v2t(z);

    % calculate error
    e = t2v(inv(zt) * (inv(xt1) * xt2));

    d = x2 - x1;

    % calculate intermediate jacobian
    s = sin(x1(3));
    c = cos(x1(3));
    Atmp = [-c, -s, -s * d(1) + c * d(2);
             s, -c, -c * d(1) - s * d(2)];

    % retrieve rotation matrix of z
    Rz = zt(1:2,1:2)';

    A = [Rz * Atmp;
         0, 0, -1];

    Btmp = [-Atmp(1,1:2), 0;
            -Atmp(2,1:2), 0];
    B = [Rz * Btmp;
         0, 0, 1];
end;
