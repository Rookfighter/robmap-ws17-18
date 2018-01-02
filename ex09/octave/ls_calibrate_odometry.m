% this function solves the odometry calibration problem
% given a measurement matrix Z.
% We assume that the information matrix is the identity
% for each of the measurements
% Every row of the matrix contains
% z_i = [u'x, u'y, u'theta, ux, uy, ytheta]
% Z:	The measurement matrix
% X:	the calibration matrix
% returns the correction matrix X
function X = ls_calibrate_odometry(Z)
    N = size(Z, 1);
    % max iterations
    M = 10;
    % Error is normal distributed: omega = sigma = identity
    Omega = eye(3);
    % initial solution (the identity transformation)
    X = eye(3);

    U_act = Z(:, 1:3)';
    U_est = Z(:, 4:6)';

    % error values
    e = zeros(3, N);

    for k=1:M
        % init H and b with 0
        H = zeros(9,9);
        b = zeros(1,9);

        e = U_act - X * U_est;
        J = zeros(3,9);
        tmp = sum(U_est,2)';
        J(1, 1:3) = tmp;
        J(2, 4:6) = tmp;
        J(3, 7:9) = tmp;

        % % accumulate data from measurements
        % for i=1:N
        %     ei = e(:,i);
        %     J = jacobian(i, Z);
        %
        %     b = b + ei' * Omega * J;
        %     H = H + J'  * Omega * J;
        % end
        % b = b';

        % check if error is sufficiently small
        if norm(e) <= 1e-3
            break;
        end

        % calc newton step
        pk = -inv(J) * e;
        X = X + reshape(pk, [3,3]);
    end
end

% this function computes the error of the i^th measurement in Z
% given the calibration parameters
% i:	the number of the measurement
% X:	the actual calibration parameters
% Z:	the measurement matrix, each row contains first the scan-match result
%       and then the motion reported by odometry
% e:	the error of the ith measurement
function e = error_function(i, X, Z)
    U_act = Z(i, 1:3)';
    U_est = Z(i, 4:6)';

    e =  U_act - X * U_est;
end

% derivative of the error function for the ith measurement in Z
% i:	the measurement number
% Z:	the measurement matrix
% J:	the jacobian of the ith measurement
function J = jacobian(i, Z)
    U_est = Z(i, 4:6);

    J = zeros(3,9);
    J(1, 1:3) = U_est;
    J(2, 4:6) = U_est;
    J(3, 7:9) = U_est;
end
