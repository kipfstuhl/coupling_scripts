clear all
close all
clc

% author: Luca Pegolotti on 29/11/2017

% This script computes the convergence order obtained on conforming meshes
% presented in Section 3. The Poisson problem is solved on two subdomains:
% (0,0.5)x(0,1) and (0.5,0)x(0,1). The basis functions on the interface are
% Fourier basis functions. We use P2 polynomials in the left domain and P1 
% in the right domain

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex')

% load the exact solution and forcing term
run load_exact_solution_and_f.m

% dirichlet and neumann functions. The order of the edges is the same as in
% the boundary condition flags. Here we are only using homogeneous boundary
% conditions anyway.
dir_functions = @(x) [0;0;0;0];
neu_functions = @(x) [0;0;0;0];

% we solve the problem with 20, 28, ... elements in the x direction
N = [20 28 40 56 80 114 160];

% we compute the standard H1 error in both subdomains
typerror = 'H1';

globalerror = [];
brokenerror = [];
error1 = [];
error2 = [];

for n_elements = N
    
    % Compute error of the global solution
    mesh = create_mesh(0,0,1,1,n_elements,n_elements);
    
    fespace = create_fespace(mesh,'P1',[1 1 1 1]);
    [A,b] = assembler_poisson(fespace,fun,mu,dir_functions,neu_functions);
    
    totalsol = A\b;
    globalerror = [globalerror compute_H1_error(fespace,totalsol,uex,graduex)];
    disp(['Total ', typerror ,' error (full solution) = ', num2str(globalerror(end))]);
    
    n1x = n_elements/2;
    n2x = n_elements/2;
    n1y = n_elements;
    n2y = n_elements+1;
    
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
    fespace1 = create_fespace(mesh1,'P1',bc1);
    
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
    
    % we start solving the problem with variable number of basis functions
    % defined over the interface
    n_iterations = 16;
    
    % these are the extra-diagonal matrices in the global matrix
    B1 = [];
    B2 = [];
    
    brokenerrors = [];
    errors1 = [];
    errors2 = [];
    
    for i = 1:n_iterations
        disp(['n_elements = ', num2str(n_elements) ...
              ', iteration n = ', num2str(i)]);
        
        % these vectors will contain at position i the results of b(phi_i,eta),
        % where phi_i is the ith basis function of subdomain 1 or subdomain 2
        % and eta is a particular basis function defined over the interface
        b1 = zeros(n1,1);
        b2 = zeros(n2,1);
        
        % frequency of the Fourier basis
        freq = i - 1;
        
        % at the first iteration, we just consider as basis function the
        % constant function. For each next iteration we add sin and cos function
        if (i == 1)
            % compute b1 by applying constant neumann boundary conditions on
            % the right boundary (hence, the second component of the neumann
            % functions is the only non zero component)
            b1 = apply_neumann_bc(b1,fespace1,@(x) [0;1;0;0],4);
            
            % bad practice, but here the number of iterations is small. At each
            % iteration we increase the size of B1 (and B2) by adding the newly
            % computed b(phi_i,eta)
            B1 = [B1;b1'];
            
            % compute b2 by applying constant neumann boundary conditions on
            % the left boundary
            b2 = apply_neumann_bc(b2,fespace2,@(x) [0;0;0;1],4);
            
            B2 = [B2;b2'];
        else
            % here we to the same as for the constant function, but we add the
            % terms corresponding to sin and cos integrals
            
            b1 = apply_neumann_bc(b1,fespace1,@(x) [0;sin(x(2) * pi * freq);0;0],4);
            B1 = [B1;b1'];
            
            % we put b1 to zero again, as apply_neumann_bc adds the condition
            % to the vector passed as argument
            b1 = b1*0;
            
            b1 = apply_neumann_bc(b1,fespace1,@(x) [0;cos(x(2) * pi * freq);0;0],4);
            B1 = [B1;b1'];
            
            % we do the same for the right subodmain
            b2 = apply_neumann_bc(b2,fespace2,@(x) [0;0;0;sin(x(2) * pi * freq)],4);
            B2 = [B2;b2'];
            
            b2 = b2*0;
            
            b2 = apply_neumann_bc(b2,fespace2,@(x) [0;0;0;cos(x(2) * pi * freq)],4);
            B2 = [B2;b2'];
        end
        
        % number of lagrange multipliers
        n3 = size(B1,1);
        
        % build the global matrix (note that Dirichlet boundary conditions are
        % imposed on A1 and A2 directly in the assembly, and that B1' and B2'
        % have 0 value in the rows corresponding to Dirichlet boundaries)
        A = sparse([A1 sparse(n1,n2) -B1'; sparse(n2,n1) A2 B2'; -B1 B2 sparse(n3,n3)]);
        
        % build the global right handside
        f = [rhs1;rhs2;zeros(n3,1)];
        
        % solve the linear system A u = f
        sol = A\f;
        
        % divide the solution into left and right solutions
        sol1 = sol(indices1);
        sol2 = sol(indices2);
        
        % compute broken norm of the error
        err1 = compute_H1_error(fespace1,sol1,uex,graduex);
        err2 = compute_H1_error(fespace2,sol2,uex,graduex);
        errors1 = [errors1;err1];
        errors2 = [errors2;err2];
        err = sqrt(err1^2+err2^2);
        
        disp(['Total ', typerror ' error = ', num2str(err)]);
        
        brokenerrors = [brokenerrors; err];                
    end
    error1 = [error1 errors1];
    error2 = [error2 errors2];
    brokenerror = [brokenerror brokenerrors];
end

if (exist('data','dir') == 0)
    disp('error')
    mkdir data;
end
save(['data/globalerror_P1_',typerror,'_nonconf.mat'],'globalerror');
save(['data/brokenerror_P2P1_',typerror,'_nonconf.mat'],'brokenerror');
save(['data/error_P1_',typerror,'_nonconf.mat'],'error1');
save(['data/error_P2_',typerror,'_nonconf.mat'],'error2');
