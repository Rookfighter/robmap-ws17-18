function [r] = rotm(theta)
    r = [cos(theta) -sin(theta);
         sin(theta)  cos(theta)];
end
