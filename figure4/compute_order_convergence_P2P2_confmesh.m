clear all
clc

% load the exact solution and forcing term
run load_exact_solution_poisson.m

dir_functions = @(x) [0;0;0;0];
neu_functions = @(x) [0;0;0;0];

% we solve the problem with 20, 28, ... elements in the x direction
N = [20 28 40 56 80 114 160];
nrefs = length(N);

% use 15 frequencies on the interface
nfreq = 15;

brokenerror = zeros(nfreq+1,nrefs);

count = 0;
for n_elements = N
    count = count+1;
    
    n1x = n_elements/2;
    n2x = n_elements/2;
    n1y = n_elements;
    n2y = n_elements;
    
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
       
        % in the first iteration, we just consider as basis function the
        % constant function. For each next iteration we add sin and cos
        % functions
        B1 = add_row_to_coupling_matrix_poisson(B1,fespace1,2,i);
        B2 = add_row_to_coupling_matrix_poisson(B2,fespace2,4,i);
        
        % number of lagrange multipliers
        n3 = size(B1,1);
        
        % build the global matrix
        mat = [A1 sparse(n1,n2) -B1';
               sparse(n2,n1) A2  B2';
               -B1 B2 sparse(n3,n3)];  
           
        % build the global right handside
        f = [rhs1;rhs2;zeros(n3,1)];
        
        % apply bc
        [mat,f] = apply_dirichlet_bc_global_matrix_and_rhs(mat,f,...
                    {fespace1,fespace2},dir_functions);

        % solve the linear system A u = f
        sol = mat\f;
        
        % divide the solution into left and right solutions
        sol1 = sol(indices1);
        sol2 = sol(indices2);
        
        % compute broken norm of the error
        err1 = compute_H1_error(fespace1,sol1,uex,graduex);
        err2 = compute_H1_error(fespace2,sol2,uex,graduex);
        err = sqrt(err1^2+err2^2);
                
        brokenerror(i+1,count) = err;              
    end
end

% save solution
save('data_figure4/brokenerror_conf.mat','brokenerror');
