clear all
clc

% run all scripts for figure 7

% create data folder
mkdir data_figure7

% compute convergence navier stokes with partition into 5 domains
run convergence_navier_stokes_5domains.m

% plot convergence
run plot_convergence_ns.m
