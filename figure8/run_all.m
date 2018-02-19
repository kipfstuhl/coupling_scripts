clear all
clc

% run all scripts for figure 8

% create data folder
mkdir data_figure8

% compute exact solution (attention: this could take a while)
run convergence_navier_stokes_5domains.m

% compute solution on 5 subdomains
run compute_liddriven_5domains.m

% plot streamlines
run plot_streamlines.m
