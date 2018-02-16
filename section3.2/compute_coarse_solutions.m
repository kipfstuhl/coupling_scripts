clear all
close all
clc

% author: Luca Pegolotti on 11/12/2017

% This script computes the "exact" solution for the navier stokes problem

bottom_left_corner_x = 0;
bottom_left_corner_y = 0;

L = 1;
H = 1;

nex = [128];

for n = nex
    n_elements_x = n;
    n_elements_y = n;
    
    mesh = create_mesh(bottom_left_corner_x, ...
        bottom_left_corner_y, ...
        L,H,n_elements_x,n_elements_y);
    
    bc_flags = [1 1 1 1];
    
    fespace_u = create_fespace(mesh,'P2',bc_flags);
    fespace_p = create_fespace(mesh,'P1',bc_flags);
    
    f = [0;0];
    mu = 1;
    
    U = 500;
    
    dirichlet_functions = @(x) [0 0;0 0;U 0;0 0]';
    neumann_functions = @(x) [0 0;0 0;0 0;0 0]';
    
    [A,b] = assembler_steady_navier_stokes(fespace_u,fespace_p,f,mu,dirichlet_functions,...
        neumann_functions);
    
    % solve stokes problem with u = 0 to get initial guess for newton's method
    x0 = zeros(size(fespace_u.nodes,1)*2 + size(fespace_p.nodes,1),1);
    
    % solve system with newton's method
    method.name = 'newton';
    method.f = @(u) A(u)*u-b;
    method.x0 = x0;
    method.jac = @(u) build_jac_navier_stokes(A,u,fespace_u);
    method.tol = 1e-10;
    method.maxit = 100;
    
    [sol,err,it] = solve_fluid_system(A,b,fespace_u,fespace_p,method);
    plot_fe_fluid_function(sol,'U');
    pause()
    
    display(['Size of the system is ', num2str(length(sol.u1)+length(sol.u2)+length(sol.p))]);
    nsys = length(x0);
    save(['data/coarse_sol',num2str(n),'.mat'],'sol')
    save(['data/nsys',num2str(n),'_coarse.mat'],'nsys')
end