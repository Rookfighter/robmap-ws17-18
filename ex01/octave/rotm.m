function [r] = rotm(theta)
    r = [cosd(theta) -sind(theta);
         sind(theta)  cosd(theta)];
end
