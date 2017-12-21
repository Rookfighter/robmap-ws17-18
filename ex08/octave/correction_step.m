function particles = correction_step(particles, z)

% Weight the particles according to the current map of the particle
% and the landmark observations z.
% z: struct array containing the landmark observations.
% Each observation z(j) has an id z(j).id, a range z(j).range, and a bearing z(j).bearing
% The vector observedLandmarks indicates which landmarks have been observed
% at some point by the robot.

% Number of particles
numParticles = length(particles);

% Number of measurements in this time step
m = size(z, 2);

% sensor noise matrix
Q_t = [0.1 0;
       0   0.1];

% default importance weight for particle
p0 = 1 / numParticles;

% process each particle
for i = 1:numParticles
    % abbrev for estimated robot pose of this particle
    rpose = particles(i).pose;

    % process each measurement
    for j = 1:m
        % Get the id of the landmark corresponding to the j-th observation
        % particles(i).landmarks(l) is the EKF for this landmark
        l = z(j).id;
        % vector representation of the measurement
        zvec = [z(j).range; z(j).bearing];

        % The (2x2) EKF of the landmark is given by
        % its mean particles(i).landmarks(l).mu
        % and by its covariance particles(i).landmarks(l).sigma

        % If the landmark is observed for the first time:
        if particles(i).landmarks(l).observed == false

            % calc mean by inverse measurement model
            mu_j = rpose(1:2) + [zvec(1) * cos(rpose(3) + zvec(2));
                                 zvec(1) * sin(rpose(3) + zvec(2))];
            particles(i).landmarks(l).mu = mu_j;

            % get the Jacobian with respect to the landmark position
            [h, H] = measurement_model(particles(i), z(j));

            % calc covariance
            H_inv = inv(H);
            sigma_j =  H_inv * Q_t * H_inv';
            particles(i).landmarks(l).sigma = sigma_j;

            % Indicate that this landmark has been observed
            particles(i).landmarks(l).observed = true;
            % set default importance weight for particle
            particles(i).weight = p0;
        else

            % some abbrevs
            mu_j    = particles(i).landmarks(l).mu;
            sigma_j = particles(i).landmarks(l).sigma;

            % get the expected measurement
            [zexp, H] = measurement_model(particles(i), z(j));

            % calc measurement covariance
            Q     = H * sigma_j * H' + Q_t;
            Q_inv = inv(Q);

            % calc Kalman gain
            K = sigma_j * H' * Q_inv;

            % calc error between the z and expectedZ; normalize angle
            zerr    = zvec - zexp;
            zerr(2) = normalize_angle(zerr(2));

            % update the mean and covariance of the EKF for this landmark
            particles(i).landmarks(l).mu    = mu_j + K * zerr;
            particles(i).landmarks(l).sigma = (eye(2) - K * H) * sigma_j;

            % compute the likelihood of this observation, multiply with the former weight
            % to account for observing several features in one time step
            particles(i).weight = inv(sqrt(norm1(2 * pi * Q))) * exp(-0.5 * zerr' * Q_inv * zerr);
            % same but less efficient (save one inversion and squareroot):
            % sig = Q;
            % mu  = zexp;
            % w   = mvnpdf(zvec, mu, sig);
        end
    end % measurement loop
end % particle loop
end
