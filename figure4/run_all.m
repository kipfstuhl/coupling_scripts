clear all
clc

% run all scripts for figure 4

% create data folder
mkdir data_figure4

% compute global error for reference
run compute_global_error.m

% compute convergence with conforming mesh
run compute_order_convergence_P2P2_confmesh.m

% compute convergence with non-conforming mesh
run compute_order_convergence_P2P2_nonconfmesh.m

% compute convergence with non-conforming mesh and 2 gauss points
run compute_order_convergence_P2P2_nonconfmesh_2gausspoints.m

% plot convergence
run plot_convergence.m

% plot convergence
run plot_instability_gausspoints.m