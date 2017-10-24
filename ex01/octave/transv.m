function [y] = transv(x, t)
    y = rotm(t(3)) * (x + t(1:2));
end
