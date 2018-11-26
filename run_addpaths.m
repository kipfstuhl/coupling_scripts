clear all
clc

% Runs the addpaths script in feamat. This allows to use the functions in
% the library from every directory in the system tree
run feamat/addpaths.m

% Add common folder to path
addpath(genpath('./common/'))

% set this variable to the location in which GeoPDEs has been decompressed
geopdes_location = '';

% Add geoPDEs folders
addpath(genpath([geopdes_location,'/geopdes']));
addpath(genpath([geopdes_location,'/geopdes_hierarchical']));
addpath(genpath([geopdes_location,'/nurbs']));