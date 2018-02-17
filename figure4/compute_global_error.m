clear all
close all
clc

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex')

% load the exact solution and forcing term
run load_exact_solution_poisson.m

% dirichlet and neumann functions. The order of the edges is the same as in
% the boundary condition flags. Here we are only using homogeneous boundary
% conditions anyway.
dir_functions = @(x) [0;0;0;0];
neu_functions = @(x) [0;0;0;0];

% we solve the problem with 20, 28, ... elements in the x direction
N = [20 28 40 56 80 114 160];
nrefs = length(N);

nfreq = 15;

% we compute the standard H1 error in both subdomains
globalerror = zeros(nrefs,1);

count = 0;
for n_elements = N
    count = count + 1;
    % Compute error of the global solution
    mesh = create_mesh(0,0,1,1,n_elements,n_elements);
    
    fespace = create_fespace(mesh,'P2',[1 1 1 1]);
    [A,b] = assembler_poisson(fespace,fun,mu,dir_functions,neu_functions);
    
    totalsol = A\b;
    globalerror(count) = compute_H1_error(fespace,totalsol,uex,graduex);
    disp(['Total global H1 error = ', num2str(globalerror(count))]);
end

save('data_figure3/globalerror.mat','globalerror');
