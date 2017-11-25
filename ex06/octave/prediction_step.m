function [mu, sigma, sigma_points] = prediction_step(mu, sigma, u)
% Updates the belief concerning the robot pose according to the motion model.
% mu: state vector containing robot pose and poses of landmarks obeserved so far
% Current robot pose = mu(1:3)
% Note that the landmark poses in mu are stacked in the order by which they were observed
% sigma: the covariance matrix of the system.
% u: odometry reading (r1, t, r2)
% Use u.r1, u.t, and u.r2 to access the rotation and translation values

% For computing lambda.
global scale;

% Compute sigma points
sigma_points = compute_sigma_points(mu, sigma);

% Dimensionality
n = length(mu);
% lambda
lambda = scale - n;

% transform matrix
% general remark: it would be more efficient to index mu(1:3)
% and sigma(1:3,1:3) than using the F matrix. map is not affected
% in this step
F = zeros(3, n);
F(1:3,1:3) = eye(3);

% define motion model
g = @(xk) xk + F' * [u.t * cos(xk(3) + u.r1);
                     u.t * sin(xk(3) + u.r1);
                     u.r1 + u.r2 ];

% calc transformed sigma points
% TODO any idea how to vectorize this?
y = zeros(size(sigma_points));
for i = 1:length(sigma_points)
    y(:,i) = g(sigma_points(:,i));
end
% normalize angles of transformed sigma points
y(3,:) = normalize_angle(y(3,:));

% Computing the weights for recovering the mean
wm = [lambda/scale, repmat(1/(2*scale),1,2*n)];
wc = wm;

% calc estimated mu
% TODO why use cos and sin? => same result as normal average!
mu = y * wm';
mu(3) = normalize_angle(mu(3));

% trig = [sin(y(3,:));
%         cos(y(3,:))];
% trig = trig * wm';
% mu(3,:) = atan2(trig(1), trig(2));

% replicate mu for vectorized multiplication
mu_r = repmat(mu, 1, length(sigma_points));
% calc difference of sigma points and mu
diff = sigma_points - mu_r;
diff(3,:) = normalize_angle(diff(3,:));
% calc estimated sigma
sigma  = diff * diag(wc) * diff';

% Motion noise
motionNoise = 0.1;
R3 = [motionNoise, 0, 0;
     0, motionNoise, 0;
     0, 0, motionNoise/10];
R = zeros(size(sigma,1));
R(1:3,1:3) = R3;

% add motion noise to sigma
sigma = sigma + R;

end
