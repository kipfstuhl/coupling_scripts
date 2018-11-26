clear all
clc

% run all scripts for figure 10

% create data folder
mkdir data_figure10

run plot_solution.m

% in order to compute the convergence orders wrt to the fine solution, we
% first interpolate the fine solution on all the coarse meshes (to speed up
% the computation of the error). This is not completely correct from a
% theoretical point of view but the approximation does not prevent to reach
% the desired convergence rates
interpolate_on_coarse_spaces_bifurcation('data_figure10',1:6);

% compute errors for convergence
run compute_order_convergence_FEM_FEM.m

% plot convergence rates
run plot_convergence.m