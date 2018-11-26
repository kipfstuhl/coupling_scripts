clear all
clc

% run all scripts for figure 13

% create folder for data
mkdir data_figure14

% interpolate fine finite element 
interpolate_on_coarse_spaces_bifurcation('data_figure14',2,'_fem_iga')

% compute convergences
run compute_convergence_fem_iga_P3P2.m

% plot convergence in Omega_1 and Omega_2
run plot_convergence.m


