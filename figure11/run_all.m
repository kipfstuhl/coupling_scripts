clear all
clc

% run all scripts for figure 11

% create data folder
mkdir data_figure11

% we interpolate the exact solution on the coarse meshes
interpolate_on_coarse_spaces_bifurcation('data_figure11',[2 6]);

% plot errors for each refinement level and number of frequencies
for ref = [2 6]
    for freq = [0 2 4]
        figure
        plot_error_bifurcation(ref,freq,3,-18.4207,-0.6);
    end
end