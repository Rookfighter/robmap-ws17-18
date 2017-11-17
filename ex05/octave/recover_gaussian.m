function [mu, sigma] = recover_gaussian(sigma_points, w_m, w_c)
% This function computes the recovered Gaussian distribution (mu and sigma)
% given the sigma points (size: nx2n+1) and their weights w_m and w_c:
% w_m = [w_m_0, ..., w_m_2n], w_c = [w_c_0, ..., w_c_2n].
% The weight vectors are each 1x2n+1 in size,
% where n is the dimensionality of the distribution.

% compute new mu
mu = sigma_points * w_m';

% some help variables
n = length(mu);
mu_r = repmat(mu, 1, 2*n+1);
% define intermediate result
tmp = sigma_points - mu_r;

% solve iterative
% sigma = zeros(n, n);
% for i = 1:2*n+1
%     d = sigma_points(:,i) - mu;
%     sigma += w_c(i) * d *  d';
% end

% compute new sigma
sigma  = tmp * diag(w_c) * tmp';

end
