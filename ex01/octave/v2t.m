% v2t.m
%
%     Author: Fabian Meyer
% Created On: 21 Oct 2017

function [t] = v2t(x)
    t = [cos(x(3)) -sin(x(3)) x(1);
         sin(x(3))  cos(x(3)) x(2);
         0          0         1];
end
