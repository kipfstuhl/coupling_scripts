clear all
close all
clc

mesh = read_mesh('../meshes/bifurcation_smooth_distorted_fine.msh');
disp('Mesh has been read');

fespace_u = create_fespace(mesh,'P3',[1 0 1 0 1 1]);
fespace_p = create_fespace(mesh,'P2',[]);

n_nodes_u = size(fespace_u.nodes,1);

disp(['# dofs velocity:',num2str(n_nodes_u)]);

r_in = 0.5;
U = 1;

dir_in = @(x) [0 0; 0 0; 0 0;0 0; 0 0; (r_in^2 - x(2)^2)/r_in^2 * U 0]';
neu0 = @(x) zeros(2,6);

[A,b] = assembler_steady_navier_stokes(fespace_u,fespace_p,[0;0],1,dir_in,neu0);
% [A,b] = assembler_steady_stokes(fespace_u,fespace_p,[0;0],1,dir_in,neu0);

J = @(u) build_jac_navier_stokes(A,u,fespace_u);

u0 = zeros(length(b),1);
x0 = A(u0)\b;

% solve system with newton's method
f = @(u) A(u)*u-b;
tol = 1e-8;
maxit = 20;

method.name = 'newton';
method.f = f;
method.jac = J;
method.tol = tol;
method.maxit = maxit;
method.x0 = x0;
[sol,err,it] = solve_fluid_system(A,b,fespace_u,fespace_p,method);
% [sol,err,it] = solve_fluid_system(A,b,fespace_u,fespace_p);

save('solutions/exsol.mat','sol');
export_vtk_fluid(sol,'fine_sol');