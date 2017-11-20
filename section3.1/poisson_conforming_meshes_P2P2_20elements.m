clear all
clc

% author: Luca Pegolotti on 20/11/2017

% This script performs the numerical simulations regarding conforming meshes 
% in Section 3.1. The Poisson problem is solved on two subdomains: (0,0.5)x(0,1) 
% and (0.5,0)x(0,1). The basis functions on the interface are Fourier basis
% functions.

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex')

% define the exact solution
alpha = 100;

uex = @(x) alpha * (1 - x(1,:)) .* x(1,:).^(1/1) .* (1 - x(2,:)) .* x(2,:) .* sin(1/3 - x(1,:) .* x(2,:).^2);

uexdx = @(x,y) - alpha * x .* y .* sin(x .* y.^2 - 1/3) .* (y - 1) - ...
                 alpha * y .* sin(x .* y.^2 - 1/3) .* (x - 1) .* (y - 1) - ...
                 alpha * x .* y.^3 .* cos(x .* y.^2 - 1/3) .* (x - 1) .* (y - 1);
uexdy = @(x,y) - alpha * x * y * sin(x*y^2 - 1/3) * (x - 1) - ... 
                 alpha * x * sin(x * y^2 - 1/3) * (x - 1) * (y - 1) - ...
                 2 * alpha * x^2 * y^2 * cos(x * y^2 - 1/3) * (x - 1) * (y - 1);
graduex = @(x) [uexdx(x(1,:),x(2,:));uexdy(x(1,:),x(2,:))];

% define forcing term
funxy = @(x,y) 2 * alpha * x * sin(x * y^2 - 1/3) * (x - 1) + ... 
               2 * alpha * y * sin(x * y^2 - 1/3) * (y - 1) + ...
               2 * alpha * x * y^3 * cos(x * y^2 - 1/3) * (y - 1) + ...
               4 * alpha * x^2 * y^2 * cos(x * y^2 - 1/3) * (x - 1) + ...
               2 * alpha * y^3 * cos(x * y^2 - 1/3) * (x - 1) * (y - 1) + ...
               6 * alpha * x^2 * y * cos(x * y^2 - 1/3) * (x - 1) * (y - 1) - ...
               alpha * x * y^5 * sin(x*y^2 - 1/3) * (x - 1)*(y - 1) - ...
               4 * alpha * x^3 * y^3 * sin(x * y^2 - 1/3) * (x - 1) * (y - 1);
fun = @(x) funxy(x(1),x(2));

% define diffusivity
mu_const = 1;
mu = @(x) mu_const;

% discretization parameters: n_elements in y direction for both
% subdomains, n_elements/2 elements in x direction for both subdomains
n_elements = 20;
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

% dirichlet and neumann functions. The order of the edges is the same as in 
% the boundary condition flags. Here we are only using homogeneous boundary 
% conditions anyway.
dir_functions = @(x) [0;0;0;0];
neu_functions = @(x) [0;0;0;0];

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
n_iterations = 3;

% these are the extra-diagonal matrices in the global matrix
B1 = [];
B2 = [];

% before the computation, we create fine finite element spaces that are
% used only for the visualization (since the function plot_on_fespace
% uses only the degrees of freedom of the vertices of triangles, we first
% interpolate the quadratic solution on a fine mesh and visualize that
% instead)
n_elements_visualization = 50;
finemesh1 = create_mesh(xp1,yp1,L1,H1,n_elements_visualization/2,n_elements_visualization);
finefespace1 = create_fespace(finemesh1,'P2',[1 1 1 1]);
finemesh2 = create_mesh(xp2,yp2,L2,H2,n_elements_visualization/2,n_elements_visualization);
finefespace2 = create_fespace(finemesh2,'P2',[1 1 1 1]);

for i = 1:n_iterations
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
        b1 = apply_neumann_bc(fespace1,b1,@(x) [0;1;0;0]);
        
        % bad practice, but here the number of iterations is small. At each
        % iteration we increase the size of B1 (and B2) by adding the newly
        % computed b(phi_i,eta)
        B1 = [B1;b1'];
        
        % compute b2 by applying constant neumann boundary conditions on
        % the left boundary 
        b2 = apply_neumann_bc(fespace2,b2,@(x) [0;0;0;1]);
        
        B2 = [B2;b2'];
    else
        % here we to the same as for the constant function, but we add the
        % terms corresponding to sin and cos integrals
        
        b1 = apply_neumann_bc(fespace1,b1,@(x) [0;sin(x(2) * pi * freq);0;0]);
        B1 = [B1;b1'];
        
        % we put b1 to zero again, as apply_neumann_bc adds the condition
        % to the vector passed as argument
        b1 = b1*0; 
        
        b1 = apply_neumann_bc(fespace1,b1,@(x) [0;cos(x(2) * pi * freq);0;0]);
        B1 = [B1;b1'];
    
        % we do the same for the right subodmain
        b2 = apply_neumann_bc(fespace2,b2,@(x) [0;0;0;sin(x(2) * pi * freq)]);
        B2 = [B2;b2'];
        
        b2 = b2*0; 
        
        b2 = apply_neumann_bc(fespace2,b2,@(x) [0;0;0;cos(x(2) * pi * freq)]);
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
    
    % create contour plot
    
    % interpolate on fine finite element spaces
    interpsol1 = interp_on_fespace(fespace1,sol1,finefespace1);
    interpsol2 = interp_on_fespace(fespace2,sol2,finefespace2);

    % draw 20 contour levels between the minimum and the maximum values at
    % the degrees of freedom
    mm = -0.55;
    MM = 1.66;
    levels = linspace(mm,MM,20);
    
    % plot contours
    figure(i)
    box
    hold on        
    set(gcf, 'renderer','painters')
    [h1]=plot_solution_on_fespace(finefespace1,interpsol1,'contourf',levels);
    [h2]=plot_solution_on_fespace(finefespace2,interpsol2,'contourf',levels);
    hold off
    
    % set axis and colours
    colorbar()
    caxis([mm MM]);
    axis square
    xlabel('$x$');
    ylabel('$y$');
    set(gca, 'fontsize',15);  
    
    % set title
    title(['$N_\Gamma$ = ', num2str(1+2*(i-1))]);
    
    % plot the derivatives at the interface and the lagrange multiplier
    figure(i+n_iterations)
    box
    hold on

    % we compute the derivatives at x = 0.5
    x0 = 0.5;
    
    % epsilon used for the numerical approximation via finite differences
    % of the derivatives
    epsil = 1e-9;

    % approximate derivative on the right subdomain
    [yy,u1] = get_values_over_line(fespace1,sol1,100,x0-epsil,'Ypar');
    [~ ,u2] = get_values_over_line(fespace1,sol1,100,x0,'Ypar');

    plot(yy,(u2-u1)/epsil,'r')
    
    % approximate derivative on the right subdomain
    [yy,u1] = get_values_over_line(fespace2,sol2,100,x0,'Ypar');
    [~ ,u2] = get_values_over_line(fespace2,sol2,100,x0+epsil,'Ypar');

    plot(yy,(u2-u1)/epsil,'b')

    % plot the exact derivative 
    plot(yy,uexdx(x0*yy.^0,yy),'k')
    
    % plot the Lagrange multiplier
    coeffs = sol(n1+n2+1:end);
    urec = yy*0 + coeffs(1);
    coeffs = coeffs(2:end);

    for j = 1:freq
        urec = urec + coeffs(2*j-1)*sin(pi*j*yy)+coeffs(2*j)*cos(pi*j*yy);
    end
    plot(yy,urec/mu_const,'color',[0 147 23]/255)

    % set axis and colours
    axis square
    axis([0 1 -3 1])
    xlabel('$y$');
    ylabel('$\partial y/\partial n$');
    set(gca, 'fontsize',15); 
    
    % set legend
    legend('Domain 1', 'Domain 2', 'Exact solution', 'Lagrange multiplier');
    
    % set title
    title(['$N_\Gamma$ = ', num2str(1+2*(i-1))]);
        
end



