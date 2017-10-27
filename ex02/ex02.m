% ex02.m
%
%     Author: Fabian Meyer
% Created On: 27 Oct 2017

clc;
clear all;
close all;

% constants for indexing
COLORED = 1;
BLANK = 2;

% p(x_t+1 | x_t, u_t)
% state transition distribution
p_std = [1 0.9;
         0 0.1];

% p(z_t | x_t)
% measurement distribution
p_md = [0.7 0.2;
        0.3 0.8];

% initial believe is uniform distributed
% robot does not know if object is colored or not
bel_0 = [0.5; 0.5];

% prediction of bayes filter
pred_1 = [0; 0];
pred_1(BLANK) = p_std(BLANK, BLANK) * bel_0(BLANK) + p_std(BLANK, COLORED) * bel_0(COLORED);
pred_1(COLORED) = p_std(COLORED, BLANK) * bel_0(BLANK) + p_std(COLORED, COLORED) * bel_0(COLORED);

% normalizer
n_1 = 1 / (p_md(COLORED, BLANK) * pred_1(BLANK) + p_md(COLORED, COLORED) * pred_1(COLORED));
% final belief of bayes filter
bel_1 = [0; 0];
bel_1(BLANK) = n_1 * p_md(COLORED, BLANK) * pred_1(BLANK);
bel_1(COLORED) = n_1 * p_md(COLORED, COLORED) * pred_1(COLORED);

bel_0
pred_1
n_1
bel_1
