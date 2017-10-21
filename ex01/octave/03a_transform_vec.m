% hcoord_test.m
%
%     Author: Fabian Meyer
% Created On: 21 Oct 2017

clc;
clear all;
close all;

addpath('tools');

x1 = [1.0; 2.0; pi];
x2 = [2.0; 3.0; 0.5*pi];
x3 = [0.5; 1.5; 0.25*pi];
x4 = [0.5; 1.5; 0.25*pi];
xtot = x1 + x2 + x3 + x4;
xtot(3) = normalize_angle(xtot(3));
xrot = [cos(xtot(3)) -sin(xtot(3));
        sin(xtot(3))  cos(xtot(3))];

t1 = v2t(x1);
t2 = v2t(x2);
t3 = v2t(x3);
t4 = v2t(x4);
ttot = t1 * t2 * t3 * t4;
xt = t2v(ttot);

xtot
xrot
xt
ttot