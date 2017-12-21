function [n] = norm1(A)
% Calculates L1 norm of A.

tmp = abs(A);
n = sum(sum(tmp));

end
