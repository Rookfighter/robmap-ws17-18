% computes a calibrated vector of odometry measurements
% by applying the bias term to each line of the measurements
% X: 	3x3 matrix obtained by the calibration process
% U: 	Nx3 matrix containing the odometry measurements
% C:	Nx3 matrix containing the corrected odometry measurements

function C = apply_odometry_correction(X, U)
    X = [X(1:3)';
         X(4:6)';
         X(7:9)'];

    C = U * X';
end
