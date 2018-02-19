clear all
clc

% run all scripts for figure 5

% create data folder
mkdir data_figure5

% compute approximation of inf-sup constant
run compute_beta.m

% plot beta
run plot_beta.m