% 03b_landmark.m
%
%     Author: Fabian Meyer
% Created On: 21 Oct 2017

clc;
clear all;
close all;

xt = [1; 1; pi / 2];
z = [2; 0];

% transform rel landmark and robot pose to homogenous coords
zhc = [z; 1];
xthc = v2t(xt);

% calc absolute landmark pos
lmhc = xthc * zhc;
lm = [lmhc(1) / lmhc(3) lmhc(2) / lmhc(3)] ;

lmhc
lm
