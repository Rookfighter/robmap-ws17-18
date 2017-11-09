function [mu, sigma] = prediction_step(mu, sigma, u)
% Updates the belief concerning the robot pose according to the motion model,
% mu: 2N+3 x 1 vector representing the state mean
% sigma: 2N+3 x 2N+3 covariance matrix
% u: odometry reading (r1, t, r2)
% Use u.r1, u.t, and u.r2 to access the rotation and translation values

% define some abbrevs
n = size(mu, 1);

% define f matrix as sparse
F = sparse(3, n);
F(1:3,1:3) = eye(3);

% define motion noise
motionNoise = 0.1;
R3 = [motionNoise, 0, 0; 
      0, motionNoise, 0; 
      0, 0, motionNoise/10];
R = sparse(n);
R(1:3,1:3) = R3;

% calc odometry motion model
motion = [u.t * cos(mu(3) + u.r1);
          u.t * sin(mu(3) + u.r1);
          u.r1 + u.r2 ];

% calc jacobian for prediction step
Gx = [0 0 -u.t * sin(mu(3) + u.r1);
      0 0  u.t * cos(mu(3) + u.r1);
      0 0  0];
G = speye(n) + F' * Gx * F;

% calc estimated mu and sigma 
mu = mu + F' * motion;
sigma = G * sigma * G' + F' * R * F;

mu(3) = normalize_angle(mu(3));

end
