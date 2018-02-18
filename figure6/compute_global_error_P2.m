clear all
clc

% load the exact solution and forcing term
run load_exact_solution_poisson.m

% dirichlet and neumann functions
dir_functions = @(x) [0;0;0;0];
neu_functions = @(x) [0;0;0;0];

% we solve the problem with 20, 28, ... elements in the x direction
N = [20 28 40 56 80 114 160];
nrefs = length(N);

nfreq = 15;

globalerror = zeros(nrefs,1);

count = 0;
for n_elements = N
    count = count + 1;
    
    % create mesh and finite element space
    mesh = create_mesh(0,0,1,1,n_elements,n_elements);
    fespace = create_fespace(mesh,'P2',[1 1 1 1]);
    
    % assemble matrix and right hand side
    [A,b] = assembler_poisson(fespace,fun,mu,dir_functions,neu_functions);
    
    % solve problem
    totalsol = A\b;
    globalerror(count) = compute_H1_error(fespace,totalsol,uex,graduex);
    disp(['Total global H1 error = ', num2str(globalerror(count))]);
end

% save solution
save('data_figure5/globalerrorP2.mat','globalerror');
