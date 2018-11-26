clear all
clc

% Runs the addpaths script in feamat. This allows to use the functions in
% the library from every directory in the system tree
run feamat/addpaths.m

% Add common folder to path
addpath(genpath('./common/'))
addpath(genpath('./nvr_stokes_iga'))

