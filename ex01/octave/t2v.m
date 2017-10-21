% t2v.m
%
%     Author: Fabian Meyer
% Created On: 21 Oct 2017

function [x] = t2v(t)
    x = [t(1,3) / t(3,3);
         t(2,3) / t(3,3);
         acos(t(1,1) / t(3,3))];
end
