clear all
clc

% run all scripts for figure 13

% create folder for data
mkdir data_figure13

% run simulations when using P2P1 in domain 3
run compute_convergence_fem_iga_P2P1.m

% run simulations when using P3P2 in domain 3
run compute_convergence_fem_iga_P3P2.m

% plot convergence
run plot_convergence.m
