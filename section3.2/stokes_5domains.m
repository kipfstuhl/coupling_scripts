clear all
close all
clc

% author: Luca Pegolotti on 28/11/2017

% This script performs the numerical simulations on the Stokes problem (lid-cavity problem)
% with non-conforming meshes and a partition of Omega into 3 subdomains

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex')
 
h_fine = 1/100;

% define parameters and boundary conditions
U = 10;
f = [0;0];
nu = 1;
dirichlet_functions = @(x) [0 0;0 0;U*(x(2) == 1) 0;0 0]';
neumann_functions = @(x) [0 0;0 0;0 0;0 0]';

% First solve on fine mesh to get a reference solution
bottom_left_corner_x = 0;
bottom_left_corner_y = 0;

L = 1;
H = 1;

n_elements_x = 120;
n_elements_y = 120;

mesh = create_mesh(0,0,L,H,n_elements_x,n_elements_y);

bc_flags = [1 1 1 1];

fespace_u = create_fespace(mesh,'P2',bc_flags);
fespace_p = create_fespace(mesh,'P1',bc_flags);

[A,b] = assembler_steady_stokes(fespace_u,fespace_p,f,nu,dirichlet_functions,neumann_functions);

refsol = solve_fluid_system(A,b,fespace_u,fespace_p);

% create the mesh and fespaces for domain 1
xp1 = 0;
yp1 = 0.3;
L1 = 1;
H1 = 0.7;

n1x = 100;
n1y = n1x*H1;
mesh1 = create_mesh(xp1,yp1,L1,H1,n1x,n1y);

h_coarse = L1/n1x;

bc_flags = [0 1 1 1];
fespace1_u = create_fespace(mesh1,'P2',bc_flags);
fespace1_p = create_fespace(mesh1,'P1',bc_flags);

% create the mesh and fespaces for domain 2
xp2 = 0.7;
yp2 = 0;
L2 = 0.3;
H2 = 0.3;

n2x = floor(L2/h_fine);
n2y = floor(H2/h_fine);
mesh2 = create_mesh(xp2,yp2,L2,H2,n2x,n2y);

bc_flags = [1 1 0 0];
fespace2_u = create_fespace(mesh2,'P2',bc_flags);
fespace2_p = create_fespace(mesh2,'P1',bc_flags);

% create the mesh and fespaces for domain 3
xp3 = 0.15;
yp3 = 0;
L3 = 0.55;
H3 = 0.3;

n3x = floor(L3/h_coarse);
n3y = floor(H3/h_coarse);
mesh3 = create_mesh(xp3,yp3,L3,H3,n3x,n3y);

bc_flags = [1 0 0 0];
fespace3_u = create_fespace(mesh3,'P2',bc_flags);
fespace3_p = create_fespace(mesh3,'P1',bc_flags);

% create the mesh and fespaces for domain 4
xp4 = 0;
yp4 = 0;
L4 = 0.15;
H4 = 0.15;

n4x = floor(L4/h_fine);
n4y = floor(H4/h_fine);
mesh4 = create_mesh(xp4,yp2,L4,H4,n4x,n4y);

bc_flags = [1 0 0 1];
fespace4_u = create_fespace(mesh4,'P2',bc_flags);
fespace4_p = create_fespace(mesh4,'P1',bc_flags);

% create the mesh and fespaces for domain 5
xp5 = 0;
yp5 = 0.15;
L5 = 0.15;
H5 = 0.15;

n5x = n4x;
n5y = n4y;
mesh5 = create_mesh(xp5,yp5,L5,H5,n5x,n5y);

bc_flags = [0 0 0 1];
fespace5_u = create_fespace(mesh5,'P2',bc_flags);
fespace5_p = create_fespace(mesh5,'P1',bc_flags);
%%
meshes = {};
meshes{end+1} = mesh1;
meshes{end+1} = mesh2;
meshes{end+1} = mesh3;
meshes{end+1} = mesh4;
meshes{end+1} = mesh5;

draw_multimesh(meshes)
%%

% build matrices and righ handsides for the 3 domains
[A1,rhs1] = assembler_steady_stokes(fespace1_u,fespace1_p,f,nu,dirichlet_functions,neumann_functions);
[A2,rhs2] = assembler_steady_stokes(fespace2_u,fespace2_p,f,nu,dirichlet_functions,neumann_functions);
[A3,rhs3] = assembler_steady_stokes(fespace3_u,fespace3_p,f,nu,dirichlet_functions,neumann_functions);
[A4,rhs4] = assembler_steady_stokes(fespace4_u,fespace4_p,f,nu,dirichlet_functions,neumann_functions);
[A5,rhs5] = assembler_steady_stokes(fespace5_u,fespace5_p,f,nu,dirichlet_functions,neumann_functions);

% store number of degrees of freedom for one component the velocity
n1u = size(fespace1_u.nodes,1);
n2u = size(fespace2_u.nodes,1);
n3u = size(fespace3_u.nodes,1);
n4u = size(fespace4_u.nodes,1);
n5u = size(fespace5_u.nodes,1);

% store number of degrees of freedom for the pressure
n1p = size(fespace1_p.nodes,1);
n2p = size(fespace2_p.nodes,1);
n3p = size(fespace3_p.nodes,1);
n4p = size(fespace4_p.nodes,1);
n5p = size(fespace5_p.nodes,1);


%% Build coupling blocks

n_iterations = 6;
gausspoints = 4;

% interface 1
fespaces_u = {fespace1_u, fespace2_u,fespace3_u,fespace5_u};
fespaces_p = {fespace1_p, fespace2_p,fespace3_p,fespace5_p};
bcs_flags = [1 0 0 0; 0 0 1 0; 0 0 1 0; 0 0 1 0]';
base_freq = 1/L1;
label = 'xpar';
n_iterations1 = n_iterations;
[blocks_interface1,blocks_interface1_t] = couple_navier_stokes_solutions(fespaces_u,fespaces_p,bcs_flags,base_freq,n_iterations1,label,gausspoints);
nlagmul1 = (n_iterations1-1)*4 + 2;

% interface 2
fespaces_u = {fespace4_u, fespace5_u};
fespaces_p = {fespace4_p, fespace5_p};
bcs_flags = [0 0 1 0; 1 0 0 0]';
base_freq = 1/L4;
label = 'xpar';
n_iterations2 = n_iterations;
[blocks_interface2,blocks_interface2_t] = couple_navier_stokes_solutions(fespaces_u,fespaces_p,bcs_flags,base_freq,n_iterations2,label,gausspoints);
nlagmul2 = (n_iterations2-1)*4 + 2;

% interface 3
fespaces_u = {fespace3_u, fespace4_u, fespace5_u};
fespaces_p = {fespace3_p, fespace4_p, fespace5_p};
bcs_flags = [0 0 0 1; 0 1 0 0; 0 1 0 0]';
base_freq = 1/H3;
label = 'ypar';
n_iterations3 = n_iterations;
[blocks_interface3,blocks_interface3_t] = couple_navier_stokes_solutions(fespaces_u,fespaces_p,bcs_flags,base_freq,n_iterations3,label,gausspoints);
nlagmul3 = (n_iterations3-1)*4 + 2;

% interface 4
fespaces_u = {fespace2_u, fespace3_u};
fespaces_p = {fespace2_p, fespace3_p};
bcs_flags = [0 0 0 1; 0 1 0 0]';
base_freq = 1/H3;
label = 'ypar';
n_iterations4 = n_iterations;
[blocks_interface4,blocks_interface4_t] = couple_navier_stokes_solutions(fespaces_u,fespaces_p,bcs_flags,base_freq,n_iterations4,label,gausspoints);
nlagmul4 = (n_iterations4-1)*4 + 2;


%%

n1 = 2*n1u+n1p;
n2 = 2*n2u+n2p;
n3 = 2*n3u+n3p;
n4 = 2*n4u+n4p;
n5 = 2*n5u+n5p;

N = n1+n2+n3+n4+n5;

indices1 = 1:n1;
indices2 = n1+1:n1+n2;
indices3 = n1+n2+1:n1+n2+n3;
indices4 = n1+n2+n3+1:n1+n2+n3+n4;
indices5 = n1+n2+n3+n4+1:n1+n2+n3+n4+n5;
     
A = blkdiag(A1,A2,A3,A4,A5);
     
B1 = [-blocks_interface1{1} blocks_interface1{2} blocks_interface1{3} sparse(nlagmul1,n4) blocks_interface1{4}];
B2 = [sparse(nlagmul2,n1) sparse(nlagmul2,n2) sparse(nlagmul2,n3) blocks_interface2{1} -blocks_interface2{2}];
B3 = [sparse(nlagmul3,n1) sparse(nlagmul3,n2) blocks_interface3{1} -blocks_interface3{2} -blocks_interface3{3}];
B4 = [sparse(nlagmul4,n1) blocks_interface4{1} -blocks_interface4{2} sparse(nlagmul4,n4) sparse(nlagmul4,n5)];

B = [B1;B2;B3;B4];

B1_t = [-blocks_interface1_t{1};blocks_interface1_t{2};blocks_interface1_t{3};sparse(n4,nlagmul1);blocks_interface1_t{4}];
B2_t = [sparse(n1,nlagmul2);sparse(n2,nlagmul2);sparse(n3,nlagmul2);blocks_interface2_t{1};-blocks_interface2_t{2}];
B3_t = [sparse(n1,nlagmul3);sparse(n2,nlagmul3);blocks_interface3_t{1};-blocks_interface3_t{2};-blocks_interface3_t{3}];
B4_t = [sparse(n1,nlagmul4);blocks_interface4_t{1};-blocks_interface4_t{2};sparse(n4,nlagmul4);sparse(n5,nlagmul4)];

B_t = [B1_t B2_t B3_t B4_t];

totalagmul = nlagmul1+nlagmul2+nlagmul3+nlagmul4;

zerosp = sparse(totalagmul,totalagmul);
mat = [A B_t; B zerosp];

zeroslg = zeros(totalagmul,1);
rhs = [rhs1;rhs2;rhs3;rhs4;rhs5;zeroslg];

sol = mat\rhs;

sol1 = sol(indices1);
sol2 = sol(indices2);
sol3 = sol(indices3);
sol4 = sol(indices4);
sol5 = sol(indices5);

solstr1.u1 = sol1(1:n1u);
solstr1.u2 = sol1(n1u+1:n1u*2);
solstr1.p  = sol1(2*n1u+1:end);
solstr1.fespace_u = fespace1_u;
solstr1.fespace_p = fespace1_p;

solstr2.u1 = sol2(1:n2u);
solstr2.u2 = sol2(n2u+1:n2u*2);
solstr2.p  = sol2(2*n2u+1:end);
solstr2.fespace_u = fespace2_u;
solstr2.fespace_p = fespace2_p;

solstr3.u1 = sol3(1:n3u);
solstr3.u2 = sol3(n3u+1:n3u*2);
solstr3.p  = sol3(2*n3u+1:end);
solstr3.fespace_u = fespace3_u;
solstr3.fespace_p = fespace3_p;

solstr4.u1 = sol4(1:n4u);
solstr4.u2 = sol4(n4u+1:n4u*2);
solstr4.p  = sol4(2*n4u+1:end);
solstr4.fespace_u = fespace4_u;
solstr4.fespace_p = fespace4_p;

solstr5.u1 = sol5(1:n5u);
solstr5.u2 = sol5(n5u+1:n5u*2);
solstr5.p  = sol5(2*n5u+1:end);
solstr5.fespace_u = fespace5_u;
solstr5.fespace_p = fespace5_p;
%%
plot_fe_fluid_function(solstr1,'U',[0 U]);
hold on
plot_fe_fluid_function(solstr2,'U',[0 U]);
plot_fe_fluid_function(solstr3,'U',[0 U]);
plot_fe_fluid_function(solstr4,'U',[0 U]);
hold on
plot_fe_fluid_function(solstr5,'U',[0 U]);
axis([0 1 0 1])

%% 
aux = solstr4.u1*0;

aux = apply_dirichlet_bc_rhs(aux,solstr4.fespace_u,@(x) [1;1;1;1]);

indices = find(aux == 1);

A4(indices,:)

B41 = blocks_interface2{1}';
B42 = -blocks_interface3{2}';
%% Plot stress at interface 1
yp = yp1;
epsi = 1e-8;

lms = sol(N+1:N+nlagmul1);

[x1,u11] = get_values_over_line(solstr1.fespace_u,solstr1.u2,100,yp,'Xpar');
[~ ,u21] = get_values_over_line(solstr1.fespace_u,solstr1.u2,100,yp+epsi,'Xpar');
[~ ,p1] = get_values_over_line(solstr1.fespace_p,solstr1.p,100,yp,'Xpar');

[x5,u15] = get_values_over_line(solstr5.fespace_u,solstr5.u2,100,yp,'Xpar');
[~ ,u25] = get_values_over_line(solstr5.fespace_u,solstr5.u2,100,yp-epsi,'Xpar');
[~ ,p5] = get_values_over_line(solstr5.fespace_p,solstr5.p,100,yp,'Xpar');

[x3,u13] = get_values_over_line(solstr3.fespace_u,solstr3.u2,100,yp,'Xpar');
[~ ,u23] = get_values_over_line(solstr3.fespace_u,solstr3.u2,100,yp-epsi,'Xpar');
[~ ,p3] = get_values_over_line(solstr3.fespace_p,solstr3.p,100,yp,'Xpar');

[x2,u12] = get_values_over_line(solstr2.fespace_u,solstr2.u2,100,yp,'Xpar');
[~ ,u22] = get_values_over_line(solstr2.fespace_u,solstr2.u2,100,yp-epsi,'Xpar');
[~ ,p2] = get_values_over_line(solstr2.fespace_p,solstr2.p,100,yp,'Xpar');
duy1 = (u21-u11)/epsi;
duy5 = (u15-u25)/epsi;
duy3 = (u13-u23)/epsi;
duy2 = (u12-u22)/epsi;

plot(x1,p1-duy1)
hold on
plot(x5,p5-duy5)
hold on
plot(x3,p3-duy3)
plot(x2,p2-duy2)

lagappr = lms(2)*x1.^0;

for j = 4:4:nlagmul1
    freq = j/4;
    lagappr = lagappr + lms(j)*sin(x1*pi*freq) + lms(j+2)*cos(x1*pi*freq);
    
end

plot(x1,lagappr)


% plot reference solution

% [x1,u1ref] = get_values_over_line(refsol.fespace_u,refsol.u2,100,yp-epsi/2,'Xpar');
% [~ ,u2ref] = get_values_over_line(refsol.fespace_u,refsol.u2,100,yp+epsi/2,'Xpar');
% [~ ,pref] = get_values_over_line(refsol.fespace_p,refsol.p,100,yp,'Xpar');
% 
% duref = (u2ref-u1ref)/epsi;
% plot(x1,pref-duref)

hold off


%% Plot stress at interface 2
yp = yp5;
epsi = 1e-8;

lms = sol(N+nlagmul1+1:N+nlagmul1+nlagmul2);

[x4,u14] = get_values_over_line(solstr4.fespace_u,solstr4.u2,100,yp,'Xpar');
[~ ,u24] = get_values_over_line(solstr4.fespace_u,solstr4.u2,100,yp-epsi,'Xpar');
[~ ,p4] = get_values_over_line(solstr4.fespace_p,solstr4.p,100,yp,'Xpar');

[x5,u15] = get_values_over_line(solstr5.fespace_u,solstr5.u2,100,yp,'Xpar');
[~ ,u25] = get_values_over_line(solstr5.fespace_u,solstr5.u2,100,yp+epsi,'Xpar');
[~ ,p5] = get_values_over_line(solstr5.fespace_p,solstr5.p,100,yp,'Xpar');

duy4 = (u24-u14)/epsi;
duy5 = (u15-u25)/epsi;

plot(x4,p4-duy4)
hold on
plot(x5,p5-duy5)

lagappr = lms(2)*x4.^0;

for j = 4:4:nlagmul2
    freq = j/4
    lagappr = lagappr + lms(j)*sin(x4*pi*freq/L4) + lms(j+2)*cos(x4*pi*freq/L4);
    
end

plot(x4,lagappr)

% plot reference solution

% [xref,u1ref] = get_values_over_line(refsol.fespace_u,refsol.u2,100,yp-epsi/2,'Xpar');
% [~ ,u2ref] = get_values_over_line(refsol.fespace_u,refsol.u2,100,yp+epsi/2,'Xpar');
% [~ ,pref] = get_values_over_line(refsol.fespace_p,refsol.p,100,yp,'Xpar');
% 
% duref = (u2ref-u1ref)/epsi;
% plot(xref,pref-duref+70,'--')

%% Plot stress at interface 3
xp = xp3;
epsi = 1e-8;

lms = sol(N+nlagmul1+1:N+nlagmul1+nlagmul2);

[x3,u13] = get_values_over_line(solstr3.fespace_u,solstr3.u2,100,xp,'Xpar');
[~ ,u23] = get_values_over_line(solstr3.fespace_u,solstr3.u2,100,xp+epsi,'Xpar');
[~ ,p3] = get_values_over_line(solstr3.fespace_p,solstr3.p,100,xp,'Xpar');

[x5,u15] = get_values_over_line(solstr5.fespace_u,solstr5.u2,100,yp,'Xpar');
[~ ,u25] = get_values_over_line(solstr5.fespace_u,solstr5.u2,100,yp+epsi,'Xpar');
[~ ,p5] = get_values_over_line(solstr5.fespace_p,solstr5.p,100,yp,'Xpar');

duy4 = (u24-u14)/epsi;
duy5 = (u15-u25)/epsi;

plot(x4,p4-duy4)
hold on
plot(x5,p5-duy5)

lagappr = lms(2)*x4.^0;

for j = 4:4:nlagmul2
    freq = j/4
    lagappr = lagappr + lms(j)*sin(x4*pi*freq/L4) + lms(j+2)*cos(x4*pi*freq/L4);
    
end

plot(x4,lagappr)

% plot reference solution

% [xref,u1ref] = get_values_over_line(refsol.fespace_u,refsol.u2,100,yp-epsi/2,'Xpar');
% [~ ,u2ref] = get_values_over_line(refsol.fespace_u,refsol.u2,100,yp+epsi/2,'Xpar');
% [~ ,pref] = get_values_over_line(refsol.fespace_p,refsol.p,100,yp,'Xpar');
% 
% duref = (u2ref-u1ref)/epsi;
% plot(xref,pref-duref+70,'--')


%%
export_vtk_fluid(solstr1,'sol1','U');
export_vtk_fluid(solstr2,'sol2','U');
export_vtk_fluid(solstr3,'sol3','U');
export_vtk_fluid(solstr4,'sol4','U');
export_vtk_fluid(solstr5,'sol5','U');