function [mu, sigma, observedLandmarks] = correction_step(mu, sigma, z, observedLandmarks)
% Updates the belief, i. e., mu and sigma after observing landmarks, according to the sensor model
% The employed sensor model measures the range and bearing of a landmark
% mu: 2N+3 x 1 vector representing the state mean.
% The first 3 components of mu correspond to the current estimate of the robot pose [x; y; theta]
% The current pose estimate of the landmark with id = j is: [mu(2*j+2); mu(2*j+3)]
% sigma: 2N+3 x 2N+3 is the covariance matrix
% z: struct array containing the landmark observations.
% Each observation z(i) has an id z(i).id, a range z(i).range, and a bearing z(i).bearing
% The vector observedLandmarks indicates which landmarks have been observed
% at some point by the robot.
% observedLandmarks(j) is false if the landmark with id = j has never been observed before.

% Number of measurements in this time step
m = size(z, 2);
n = size(mu, 1);

% Z: vectorized form of all measurements made in this time step: [range_1; bearing_1; range_2; bearing_2; ...; range_m; bearing_m]
% ExpectedZ: vectorized form of all expected measurements in the same form.
% They are initialized here and should be filled out in the for loop below
Z = zeros(m*2, 1);
expectedZ = zeros(m*2, 1);

% Iterate over the measurements and compute the H matrix
% (stacked Jacobian blocks of the measurement function)
% H will be 2m x 2N+3
H = zeros(2 * m, n);

% construct sensor noise matrix
Q = 0.01 * speye(2*m);

for i = 1:m
    % Get the id of the landmark corresponding to the i-th observation
	lmId = z(i).id;
    
    myz = [z(i).range;
           z(i).bearing];

    idm1 = i*2-1;
    idm2 = i*2;
    id1 = lmId*2+2;
    id2 = lmId*2+3;
    % If the landmark is obeserved for the first time:
    if (observedLandmarks(lmId)==false)

        % calc relative pos of landmark to robot
        relpos = [myz(1) * cos(myz(2) + mu(3));
                  myz(1) * sin(myz(2) + mu(3));];
        mu(id1:id2) = mu(1:2) + relpos;
        % Indicate in the observedLandmarks vector that this landmark has been observed
        observedLandmarks(lmId) = true;
    endif

    % Add the landmark measurement to the Z vector
    Z(idm1:idm2) = myz;
    
    % calc expected measurement
    delt = mu(id1:id2) - mu(1:2);
    q = delt' * delt;
    expectedZ(idm1:idm2) = [sqrt(q);
                            atan2(delt(2), delt(1)) - mu(3)];
    
    F = sparse(5,n);
    F(1:3, 1:3) = eye(3);
    F(4:5, id1:id2) = eye(2);
    
    % compute jacobian for landmark i
    Hix = (1 / q) * [ -sqrt(q) * delt(1), -sqrt(q) * delt(2),  0,  sqrt(q) * delt(1),  sqrt(q) * delt(2);
                       delt(2),           -delt(1),           -q, -delt(2),            delt(1)];
    % every Hi has two rows
    Hi = Hix * F;
    % Augment H with the new Hi
    H(idm1:idm2, :) = Hi;	
endfor

% repeat sigma and Q for suitable size
sigma_r = repmat(sigma, m);
Q_r = repmat(Q, m, 1);
% calc kalman gain
K = sigma * H' * inv(H * sigma * H' + Q);
% calc new mu and sigma
mu = mu + K * normalize_all_bearings(Z - expectedZ);
mu(3) = normalize_angle(mu(3));
sigma = (speye(n) - K * H) * sigma;

end
