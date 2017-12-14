function p = log_odds_to_prob(l)
% Convert log odds l to the corresponding probability values p.
% l could be a scalar or a matrix.

p = 1 - (1 + exp(l)).^-1;

end
