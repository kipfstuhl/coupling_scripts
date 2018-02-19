clear all
clc

% run all scripts for figure 2

% create data folder
mkdir data_figure2

% compute condition number for non-orthonormal basis
run compute_condition_number_nonortho.m

% compute condition number for orthonormal basis
run compute_condition_number_ortho.m

% plot condition numbers
run plot_cond_numbers.m