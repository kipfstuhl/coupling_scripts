clear all
close all
clc

% author: Luca Pegolotti on 28/11/2017

% This script performs the numerical simulations on the Stokes problem (lid-cavity problem)
% with non-conforming meshes and a partition of Omega into 3 subdomains

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex')



% create the mesh and fespaces for domain 1
xp1 = 0;
yp1 = 0;
L1 = 0.5;
H1 = 1;

n1x = 10;
n1y = 19;
mesh1 = create_mesh(xp1,yp1,L1,H1,n1x,n1y);

bc_flags = [1 0 1 0];
fespace1_u = create_fespace(mesh1,'P2',bc_flags);
fespace1_p = create_fespace(mesh1,'P1',bc_flags);

% create the mesh and fespaces for domain 2
xp2 = 0.5;
yp2 = 0;
L2 = 0.5;
H2 = 0.5;

n2x = 13;
n2y = 10;
mesh2 = create_mesh(xp2,yp2,L2,H2,n2x,n2y);

bc_flags = [1 1 0 0];
fespace2_u = create_fespace(mesh2,'P2',bc_flags);
fespace2_p = create_fespace(mesh2,'P1',bc_flags);

% create the mesh and fespaces for domain 3
xp3 = 0.5;
yp3 = 0.5;
L3 = 0.5;
H3 = 0.5;

n3x = 8;
n3y = 11;
mesh3 = create_mesh(xp3,yp3,L3,H3,n3x,n3y);

bc_flags = [0 1 1 0];
fespace3_u = create_fespace(mesh3,'P2',bc_flags);
fespace3_p = create_fespace(mesh3,'P1',bc_flags);

% define parameters and boundary conditions
U = 10;
f = @(x) [0;0];
nu = @(x) 1;
dirichlet_functions = @(x) [0 0;0 U*(x(1)==1);0*(x(2) == 1) 0;0 0]';

neumann_functions = @(x) [0 0;0 0;0 0;0 0]';

% build matrices and righ handsides for the 3 domains
[H1,rhs1] = assembler_steady_stokes(fespace1_u,fespace1_p,f,nu,dirichlet_functions,neumann_functions);
[H2,rhs2] = assembler_steady_stokes(fespace2_u,fespace2_p,f,nu,dirichlet_functions,neumann_functions);
[H3,rhs3] = assembler_steady_stokes(fespace3_u,fespace3_p,f,nu,dirichlet_functions,neumann_functions);

% store number of degrees of freedom for one component the velocity
n1u = size(fespace1_u.nodes,1);
n2u = size(fespace2_u.nodes,1);
n3u = size(fespace3_u.nodes,1);

% store number of degrees of freedom for one component the pressure
n1p = size(fespace1_p.nodes,1);
n2p = size(fespace2_p.nodes,1);
n3p = size(fespace3_p.nodes,1);


% build the coupling blocks

gaussorder = 7;
% blocks for the coupling over interface 1
niterations1 = 10;
nlagmult1 = (niterations1*2-1)*2;
B14 = [];
B24 = [];
B34 = [];

for i = 1:niterations1
    b1x = zeros(n1u,1);
    b2x = zeros(n2u,1);
    b3x = zeros(n3u,1);
    b1y = zeros(n1u,1);
    b2y = zeros(n2u,1);
    b3y = zeros(n3u,1);
    
    b1p = zeros(n1p,1);
    b2p = zeros(n2p,1);
    b3p = zeros(n3p,1);

    
    freq = i - 1;
    
    if (i == 1)
        b1x = apply_neumann_bc(fespace1_u,b1x,@(x) [0;1;0;0],gaussorder);
        B14 = [B14;[b1x;b1y;b1p]'];
        B14 = [B14;[b1y;b1x;b1p]'];
        
        b2x = apply_neumann_bc(fespace2_u,b2x,@(x) [0;0;0;1],gaussorder);
        B24 = [B24;[b2x;b2y;b2p]'];
        B24 = [B24;[b2y;b2x;b2p]'];
        
        b3x = apply_neumann_bc(fespace3_u,b3x,@(x) [0;0;0;1],gaussorder);
        B34 = [B34;[b3x;b3y;b3p]'];
        B34 = [B34;[b3y;b3x;b3p]'];
    else
        b1x = apply_neumann_bc(fespace1_u,b1x,@(x) [0;sin(x(2) * pi * freq);0;0],gaussorder);
        B14 = [B14;[b1x;b1y;b1p]'];
        B14 = [B14;[b1y;b1x;b1p]'];
        
        b1x = b1x*0; 
        
        b1x = apply_neumann_bc(fespace1_u,b1x,@(x) [0;cos(x(2) * pi * freq);0;0],gaussorder);
        B14 = [B14;[b1x;b1y;b1p]'];
        B14 = [B14;[b1y;b1x;b1p]'];
    
        b2x = apply_neumann_bc(fespace2_u,b2x,@(x) [0;0;0;sin(x(2) * pi * freq)],gaussorder);
        B24 = [B24;[b2x;b2y;b2p]'];
        B24 = [B24;[b2y;b2x;b2p]'];
        
        b2x = b2x*0; 
        
        b2x = apply_neumann_bc(fespace2_u,b2x,@(x) [0;0;0;cos(x(2) * pi * freq)],gaussorder);
        B24 = [B24;[b2x;b2y;b2p]'];
        B24 = [B24;[b2y;b2x;b2p]'];
        
        b3x = apply_neumann_bc(fespace3_u,b3x,@(x) [0;0;0;sin(x(2) * pi * freq)],gaussorder);
        B34 = [B34;[b3x;b3y;b3p]'];
        B34 = [B34;[b3y;b3x;b3p]'];
        
        b3x = b3x*0; 
        
        b3x = apply_neumann_bc(fespace3_u,b3x,@(x) [0;0;0;cos(x(2) * pi * freq)],gaussorder);
        B34 = [B34;[b3x;b3y;b3p]'];
        B34 = [B34;[b3y;b3x;b3p]'];
    end
end

% blocks for the coupling over interface 2
niterations2 = 5;
B25 = [];
B35 = [];
nlagmult2 = (niterations2*2-1)*2;

for i = 1:niterations2
    b2x = zeros(n2u,1);
    b3x = zeros(n3u,1);
    b2y = zeros(n2u,1);
    b3y = zeros(n3u,1);
    
    b2p = zeros(n2p,1);
    b3p = zeros(n3p,1);
    
    freq = i - 1;
    
    if (i == 1)
        b2y = apply_neumann_bc(fespace2_u,b2y,@(x) [0;0;1;0],gaussorder);
        B25 = [B25;[b2x;b2y;b2p]'];
        B25 = [B25;[b2y;b2x;b2p]'];
        
        b3y = apply_neumann_bc(fespace3_u,b3y,@(x) [1;0;0;0],gaussorder);
        B35 = [B35;[b3x;b3y;b3p]'];
        B35 = [B35;[b3y;b3x;b3p]'];
    else
        b2y = apply_neumann_bc(fespace2_u,b2y,@(x) [0;0;sin(x(1) * pi * freq/2);0],gaussorder);
        B25 = [B25;[b2x;b2y;b2p]'];
        B25 = [B25;[b2y;b2x;b2p]'];
        
        b2y = b2y*0; 
        
        b2y = apply_neumann_bc(fespace2_u,b2y,@(x) [0;0;cos(x(1) * pi * freq/2);0],gaussorder);
        B25 = [B25;[b2x;b2y;b2p]'];
        B25 = [B25;[b2y;b2x;b2p]'];
        
        b3y = apply_neumann_bc(fespace3_u,b3y,@(x) [sin(x(1) * pi * freq/2);0;0;0],gaussorder);
        B35 = [B35;[b3x;b3y;b3p]'];
        B35 = [B35;[b3y;b3x;b3p]'];
        
        b3y = b3y*0; 
        
        b3y = apply_neumann_bc(fespace3_u,b3y,@(x) [cos(x(1) * pi * freq/2);0;0;0],gaussorder);
        B35 = [B35;[b3x;b3y;b3p]'];
        B35 = [B35;[b3y;b3x;b3p]'];
    end
end

n1 = 2*n1u + n1p;
n2 = 2*n2u + n2p;
n3 = 2*n3u + n3p;


% assemble global matrix
A = sparse([H1 sparse(n1,n2) sparse(n1,n3); ...
            sparse(n2,n1) H2 sparse(n2,n3); ...
            sparse(n3,n1) sparse(n3,n2) H3]);
B = [-B14 B24 B34; sparse(nlagmult2,n1) -B25 B35];

mat = [A B'; B sparse(nlagmult1+nlagmult2,nlagmult1+nlagmult2)];
rhs = [rhs1;rhs2;rhs3;zeros(nlagmult1+nlagmult2,1)];

sol = mat\rhs;
sol1 = sol(1:n1);
sol2 = sol(n1+1:n1+n2);
sol3 = sol(n1+n2+1:n1+n2+n3);

figure(1)
plot_solution_vp(fespace1_u,fespace1_p,sol1,'U',[0 10]);
caxis([0 10])
hold on
plot_solution_vp(fespace2_u,fespace2_p,sol2,'U',[0 10]);
caxis([0 10])
hold on
plot_solution_vp(fespace3_u,fespace3_p,sol3,'U',[0 10]);
caxis([0 10])

axis([0 1 0 1])

title('Non-conforming mesh')
figure(2)

plot_solution_vp(fespace1_u,fespace1_p,sol1,'U',[0 10]);
caxis([0 10])
hold on
plot_solution_vp(fespace2_u,fespace2_p,sol2,'U',[0 10]);
caxis([0 10])
hold on
plot_solution_vp(fespace3_u,fespace3_p,sol3,'U',[0 10]);
caxis([0 10])

axis([0 1 0 1])
hold on
meshes = {}
meshes{end+1} = mesh1;
meshes{end+1} = mesh2;
meshes{end+1} = mesh3;
draw_multimesh(meshes)
hold on
title('Non-conforming mesh')

% solve full problem for comparison
% create the mesh and fespaces for domain 1
xp = 0;
yp = 0;
L = 1;
H = 1;

nx = 20;
ny = 20;
mesh = create_mesh(xp,yp,L,H,nx,ny);
fespace_u = create_fespace(mesh,'P2',[1 1 1 0]);
fespace_p = create_fespace(mesh,'P1',[0 0 0 0]);

[H,b] = assembler_steady_stokes(fespace_u,fespace_p,f,nu,dirichlet_functions,neumann_functions);
sol = H\b;
figure(3)
plot_solution_vp(fespace_u,fespace_p,sol,'U',[0 10]);
title('Full solution with conforming mesh')
