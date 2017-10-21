% 03a_transform_vec.m
%
%     Author: Fabian Meyer
% Created On: 21 Oct 2017

clc;
clear all;
close all;

addpath('tools');

% define 4 transform vecs
x1 = [1.0; 2.0; pi];
x2 = [2.0; 3.0; 0.5*pi];
x3 = [0.5; 1.5; 0.25*pi];
x4 = [0.5; 1.5; 0.25*pi];

% calc total transform vec
x = zeros(3,1);
x(1:2) = x1(1:2) + \
    rotm(x1(3)) * x2(1:2) + \
    rotm(x1(3) + x2(3)) * x3(1:2) + \
    rotm(x1(3) + x2(3) + x4(3)) * x4(1:2);
x(3) = normalize_angle(x1(3) + x2(3) + x3(3) + x4(3));

% calc homogenous transforms
t1 = v2t(x1);
t2 = v2t(x2);
t3 = v2t(x3);
t4 = v2t(x4);

% calc total homogenous transform
t = t1 * t2 * t3 * t4;
xt = t2v(t);

x
xt
t
