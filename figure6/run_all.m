clear all
clc

% run all scripts for figure 6

% create data folder
mkdir data_figure6

% compute global error with linear polynomials for reference
run compute_global_error_P1.m

% compute global error with quadratic polynomials for reference
run compute_global_error_P2.m

% compute convergence with linear polynomials in omega 1
% and quadratic polynomials in omega 2
run compute_order_convergence_P1P2_nonconfmesh.m

% plot convergence
run plot_convergence_P1P2.m