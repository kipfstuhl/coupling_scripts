clear all
clc

% run all scripts for figure 12

% options inflow (IGA)
options_inflow.element_name = 'th';     % Element type for discretization
options_inflow.degree = [ 3  3];        % Degree of the splines (pressure space)
options_inflow.regularity = [ 2  2];    % Regularity of the splines (pressure space)
options_inflow.nsub = [ 6  6];          % Number of subdivisions
options_inflow.nquad = [ 5  5];         % Points for the Gaussian quadrature rule

% options outflow 1 (FEM)
options_outflow1.polydegree_u = 'P2';
options_outflow1.polydegree_p = 'P1';
options_outflow1.ref = 1;

% options outflow 2 (IGA)
options_outflow2.element_name = 'th';
options_outflow2.degree = [ 2  2];
options_outflow2.regularity = [ 1  1];
options_outflow2.nsub = [ 5  5];
options_outflow2.nquad = [ 5  5];

[sol_fem,is] = solve_system_bifurcation_FEM_IGA(options_inflow,options_outflow1,options_outflow2,5);

% plot the solutions and meshes
plot_solution_fem_iga