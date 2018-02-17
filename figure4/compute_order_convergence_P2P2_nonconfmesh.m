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
brokenerror = zeros(nfreq+1,nrefs);

R = compute_orthonormalization_matrix(nfreq,10000);

count = 0;
for n_elements = N
    count = count+1;
    
    n1x = n_elements/2;
    n2x = n_elements/2;
    n1y = n_elements;
    n2y = n_elements + 1; % add one element to make the mesh non-conforming
    
    % create mesh on left subdomain
    xp1 = 0;
    yp1 = 0;
    L1 = 0.5;
    H1 = 1;
    
    % create a rectangular mesh with lenght L1 and height H1, with n1x elements
    % in x direction and n1y elements in y direction, with left corner being (xp1,yp1)
    mesh1 = create_mesh(xp1,yp1,L1,H1,n1x,n1y);
    
    % create mesh on right subdomain
    xp2 = L1;
    yp2 = 0;
    L2 = 0.5;
    H2 = 1;
    
    % create a rectangular mesh with lenght L2 and height H2, with n2x elements
    % in x direction and n2y elements in y direction, with left corner being (xp2,yp2)
    mesh2 = create_mesh(xp2,yp2,L2,H2,n2x,n2y);
    
    % create finite element on left subdomain. The fespace structure contains
    % information about the basis functions (in this case, polynomial), the
    % mesh and the nodes. The boundary condition flags are 1 for Dirichlet
    % boundary conditions, 0 for Neumann boundary conditions. The
    % flags are assigned in counter-clockwise order starting from the bottom
    % edge. For the left subdomain, we impose homogeneous Dirichlet at the
    % bottom, at the top and at the left boundary and we impose homogeneous
    % Neumann conditions at the right boundary.
    bc1 = [1 0 1 1];
    fespace1 = create_fespace(mesh1,'P2',bc1);
    
    % assemble matrices and rhs in the left subdomain
    [A1,rhs1] = assembler_poisson(fespace1,fun,mu,dir_functions,neu_functions);
    
    % do the same for the right subdomain (note that now we impose Neumann
    % conditions at the left boundary)
    bc2 = [1 1 1 0];
    fespace2 = create_fespace(mesh2,'P2',bc2);
    
    [A2,rhs2] = assembler_poisson(fespace2,fun,mu,dir_functions,neu_functions);
    
    % degrees of freedom in the two subdomains + indices of the degrees of
    % freedom when the dofs are stacked in a single vector
    n1 = size(A1,1);
    n2 = size(A2,1);
    
    indices1 = 1:n1;
    indices2 = n1+1:n1+n2;
    
    % these are the extra-diagonal matrices in the global matrix
    B1 = [];
    B2 = [];
        
    for i = 0:nfreq
        disp(['N elements = ',num2str(n_elements),', freq = ',num2str(i)]);
       
        % at the first iteration, we just consider as basis function the
        % constant function. For each next iteration we add sin and cos
        % function. Note that, since the mesh is non-conforming, we
        % increase the number of the gauss points for the integration from
        % 2 (default value) to 4
        B1 = add_row_to_coupling_matrix_poisson(B1,fespace1,2,i,4);
        B2 = add_row_to_coupling_matrix_poisson(B2,fespace2,4,i,4);
        
        % number of lagrange multipliers
        n3 = size(B1,1);
        
        
        % build the global matrix (note that Dirichlet boundary conditions are
        % imposed on A1 and A2 directly in the assembly, and that B1' and B2'
        % have 0 value in the rows corresponding to Dirichlet boundaries)
        mat = [A1 sparse(n1,n2) -B1';
               sparse(n2,n1) A2  B2';
               -B1 B2 sparse(n3,n3)];        
        % build the global right handside
        f = [rhs1;rhs2;zeros(n3,1)];
        
        % solve the linear system A u = f
        sol = mat\f;
        
        % divide the solution into left and right solutions
        sol1 = sol(indices1);
        sol2 = sol(indices2);
        
        % compute broken norm of the error
        err1 = compute_H1_error(fespace1,sol1,uex,graduex);
        err2 = compute_H1_error(fespace2,sol2,uex,graduex);
        err = sqrt(err1^2+err2^2);
        
        disp(['Total H1 error = ', num2str(err)]);
        brokenerror(i+1,count) = err;              
    end
end

save('data_figure3/brokenerror_nonconf.mat','brokenerror');
