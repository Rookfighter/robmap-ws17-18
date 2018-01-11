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
    maxit = 100;
    % Error is normal distributed: omega = sigma = identity
    Omega = eye(3);
    % initial solution
    X = zeros(9,1);
    % alpha param for Levenberg-Marquardt method
    % get Gauss-Newton for alpha = 0
    alpha = 1.0;


    disp('Run newton method.')
    for k=1:maxit
        % init H and b with 0
        H = zeros(9,9);
        b = zeros(1,9);

        % accumulate data from measurements
        for i=1:N
            erri = error_function(i, X, Z); %err(:,i);
            Ji   = jacobian(i, Z);

            b = b + erri' * Omega * Ji;
            H = H + Ji'  * Omega * Ji;
        end
        b = b';

        % apply Levenberg-Marquardt method
        H = H + eye(9) * alpha;
        % calc newton step
        pk = -inv(H) * b;
        X = X + pk;

        % check if we found minimum, i.e. pk is sufficiently close
        % to zero
        if norm(pk) <= 1e-3
            disp(['Newton method converged after ', num2str(k), ' iterations!'])
            break;
        end
    end
end

% this U_estfunction computes the error of the i^th measurement in Z
% given the calibration parameters
% i:	the number of the measurement
% X:	the actual calibration parameters
% Z:	the measurement matrix, each row contains first the scan-match result
%       and then the motion reported by odometry
% e:	the error of the ith measurement
function e = error_function(i, X, Z)
    U_act = Z(i, 1:3)';
    U_est = Z(i, 4:6)';

    X = [X(1:3)';
         X(4:6)';
         X(7:9)'];

    e =  U_act - X * U_est;
end

% derivative of the error function for the ith measurement in Z
% i:	the measurement number
% Z:	the measurement matrix
% J:	the jacobian of the ith measurement
function J = jacobian(i, Z)
    U_est = Z(i, 4:6);

    J = zeros(3,9);
    J(1, 1:3) = -U_est;
    J(2, 4:6) = -U_est;
    J(3, 7:9) = -U_est;
end
