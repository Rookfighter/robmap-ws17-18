% Compute the error of a pose-landmark constraint
% x 3x1 vector (x,y,theta) of the robot pose
% l 2x1 vector (x,y) of the landmark
% z 2x1 vector (x,y) of the measurement, the position of the landmark in
%   the coordinate frame of the robot given by the vector x
%
% Output
% e 2x1 error of the constraint
% A 2x3 Jacobian wrt x
% B 2x2 Jacobian wrt l
function [e, A, B] = linearize_pose_landmark_constraint(x, l, z)

    X = v2t(x);
    % retrieve rotation matrix of x
    Rx = X(1:2,1:2);

    d = l(1:2) - x(1:2);
    % calculate error
    e = Rx' * d - z;
    % jacobian of err wrt x
    s = sin(x(3));
    c = cos(x(3));
    A = [-c, -s, -s * d(1) + c * d(2);
          s, -c, -c * d(1) - s * d(2)];
    % jacobian of err wrt l
    B = -A(1:2, 1:2);

end;
