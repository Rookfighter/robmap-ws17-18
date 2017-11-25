function [mu, sigma, map] = correction_step(mu, sigma, z, map);
% Updates the belief, i.e., mu and sigma after observing landmarks,
% and augments the map with newly observed landmarks.
% The employed sensor model measures the range and bearing of a landmark
% mu: state vector containing robot pose and poses of landmarks obeserved so far.
% Current robot pose = mu(1:3)
% Note that the landmark poses in mu are stacked in the order by which they were observed
% sigma: the covariance matrix of the system.
% z: struct array containing the landmark observations.
% Each observation z(i) has an id z(i).id, a range z(i).range, and a bearing z(i).bearing
% The vector 'map' contains the ids of all landmarks observed so far by the robot in the order
% by which they were observed, NOT in ascending id order.

% For computing sigma
global scale;

% Number of measurements in this time step
m = length(z);

% Measurement noise
Q = 0.01*eye(2);

% define sensor model
h = @(x, d) [sqrt(d' * d);
             atan2(d(2), d(1)) - x(3)];

for i = 1:m

    % If the landmark is observed for the first time:
    if (isempty(find(map == z(i).id)))
        % Add new landmark to the map
        [mu, sigma, map] = add_landmark_to_map(mu, sigma, z(i), map, Q);
        % The measurement has been incorporated so we quit the correction step
        continue;
    endif

    % Compute sigma points from the predicted mean and covariance
    % This corresponds to line 6 on slide 32
    sigma_points = compute_sigma_points(mu, sigma);
    % Normalize!
    sigma_points(3,:) = normalize_angle(sigma_points(3,:));

    % Compute lambda
    n = length(mu);
    sig_len = length(sigma_points);
    lambda = scale - n;

    % extract the current location of the landmark for each sigma point
    % Use this for computing an expected measurement, i.e., applying the h function
    lm_idx = find(map==(z(i).id));
    lm_pos = [sigma_points(2*lm_idx + 2, :);
              sigma_points(2*lm_idx + 3, :)];
    rob_pos = sigma_points(1:3, :);

    % compute vector from robot to landmark
    delt = lm_pos - rob_pos(1:2, :);
    % calc expected measurements z_ex
    % TODO any idea how to vectorize this?
    z_ex = zeros(2, sig_len);
    for j = 1:sig_len
        z_ex(:, j) = h(rob_pos(:, j), delt(:, j));
    end
    % normalize bearings of z_ex
    normalize_angle(z_ex(2,:));

    % setup the weight vector for mean and covariance
    wm = [lambda/scale, repmat(1/(2*scale), 1, 2*n)];
    wc = wm;

    % calc recovered measurements
    % TODO why use cos and sin? => same result as normal average!
    trig = [sin(z_ex(2,:));
            cos(z_ex(2,:))];
    trig = trig * wm';
    zm = [z_ex(1,:) * wm';
          atan2(trig(1), trig(2))];
    % zm = z_ex * wm';
    zm(2) = normalize_angle(zm(2));

    % replicate zm for calculating diff
    zm_r = repmat(zm, 1, sig_len);
    z_diff = z_ex - zm;
    z_diff(2,:) = normalize_angle(z_diff(2,:));
    % calc innovation matrix
    S = z_diff * diag(wc) * z_diff' + Q;

    % replicate mu for calculating diff
    mu_r = repmat(mu, 1, sig_len);
    mu_diff = sigma_points - mu_r;
    mu_diff(3,:) = normalize_angle(mu_diff(3,:));
    % calc sigma_x_z
    sigma_x_z = mu_diff * diag(wc) * z_diff';

    % calc kalman gain
    K = sigma_x_z * inv(S);

    % Get the actual measurement as a vector (for computing the difference to the observation)
    z_act = [z(i).range; z(i).bearing];

    % update mu
    mu = mu + K * (z_act - zm);
    mu(3) = normalize_angle(mu(3));

    % update sigma
    sigma = sigma - K * S * K';
endfor

end
