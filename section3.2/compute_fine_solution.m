% define parameters and boundary conditions
U = 1;
f = @(x) [0;0];
nu = @(x) 1;
dirichlet_functions = @(x) [0 0;0 0;0 0;0 U*(x(1)==0)]';

neumann_functions = @(x) [0 0;0 0;0 0;0 0]';

% create the mesh and fespaces for domain 1
xp = 0;
yp = 0;
L = 1;
H = 1;

nx = 250;
ny = 250;
mesh = create_mesh(xp,yp,L,H,nx,ny);
fespace_u = create_fespace(mesh,'P2',[1 1 1 1]);
fespace_p = create_fespace(mesh,'P1',[0 0 0 0]);

[H,b] = assembler_steady_stokes(fespace_u,fespace_p,f,nu,dirichlet_functions,neumann_functions);

% manually fix value of pressure in one point of the domain to zero
nu = size(fespace_u.nodes,1);
np = size(fespace_p.nodes,1);

H(nu*2+1,:) = zeros(1,nu*2+np);
H(nu*2+1,nu*2+1) = 1;
b(nu*2+1) = 0;

sol = H\b;

save('data/full_sol.mat','sol');
save('data/fespace_u.mat','fespace_u');
save('data/fespace_p.mat','fespace_p');
%%
load('data/full_sol.mat','sol');
load('data/fespace_u.mat','fespace_u');
load('data/fespace_p.mat','fespace_p');

figure(1)
plot_solution_vp(fespace_u,fespace_p,sol,'U',[0 1]);
title('Full solution with conforming mesh')
axis([0 1 0 1])
axis square

figure(2)
plot_solution_vp(fespace_u,fespace_p,sol,'U',[0 2e-4]);
title('Full solution with conforming mesh')
axis([0.85 1 0.85 1])
axis square

