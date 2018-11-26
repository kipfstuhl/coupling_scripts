clear all
close all
clc

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex')

errsu = zeros(3,1);
errsp = zeros(3,1);

mu = 1;

dir = @(x) [exsol_u(x) exsol_u(x) exsol_u(x) exsol_u(x)  exsol_u(x)];
neu = @(x) zeros(2,6);

totalerr = [];
totalerr_lag = [];

ref_inflow = 6;
ref_out1 = 6;

% ========================= FINITE ELEMENT PART =========================

% load finite element matrices (we load the second outflow just to retrieve
% geometric information about the interface 3)
mesh_in = read_mesh(['../meshes/refinement',num2str(ref_inflow),'/inflow_distorted.msh']);
mesh_out1 = read_mesh(['../meshes/refinement',num2str(ref_out1),'/outflow1_distorted.msh']);
mesh_out2 = read_mesh(['../meshes/refinement',num2str(ref_out1),'/outflow2_distorted.msh']);
h = max(mesh_in.h,mesh_out1.h);

fespace_us = cell(2,1);
fespace_ps = cell(2,1);

% fespaces inlet
fespace_us{1} = create_fespace(mesh_in,'P2',[1 0 0 1 1]);
fespace_ps{1} = create_fespace(mesh_in,'P1',[]);

% fespaces outlet1
fespace_us{2} = create_fespace(mesh_out1,'P2',[1 0 1 0 0]);
fespace_ps{2} = create_fespace(mesh_out1,'P1',[]);

r_in = 0.5;
U = 1;
dir0 = @(x) zeros(2,6);
dir_in = @(x) [0 0; 0 0; 0 0;0 0; (r_in^2 - x(2)^2)/r_in^2 * U 0]';
neu0 = dir0;

As = cell(3,1);
bs = cell(3,1);

[As{1},bs{1}] = assembler_steady_navier_stokes(fespace_us{1},fespace_ps{1},[0;0],1,dir_in,neu0);
[As{2},bs{2}] = assembler_steady_navier_stokes(fespace_us{2},fespace_ps{2},[0;0],1,dir0,neu0);

% find interface 1
b1 = mesh_in.boundaries{2};
x1 = mesh_in.vertices(b1(end,end),1:2);
x2 = mesh_in.vertices(b1(1,1),1:2);
dif = x2 - x1;
l1 = norm(dif);
t1 = dif/l1;
n1 = [t1(2) -t1(1)];

% find interface 2
b2 = mesh_out1.boundaries{4};
x1 = mesh_out1.vertices(b2(1,1),1:2);
x2 = mesh_out1.vertices(b2(end,end),1:2);
dif = x2 - x1;
l2 = norm(dif);
t2= dif/l2;
n2 = [t2(2) -t2(1)];

% find interface 3
b3 = mesh_out2.boundaries{5};
x1 = mesh_out2.vertices(b3(1,1),1:2);
x2 = mesh_out2.vertices(b3(end,end),1:2);
dif = x2 - x1;
l3 = norm(dif);
t3= dif/l3;
n3 = [t3(2) -t3(1)];

ts = [t1;t2;t3];
ns = [n1;n2;n3];
ls = [l1 l2 l3];

nfreqs = [4 4 4];

n_nodes_us = cell(2,1);
n_nodes_ps = cell(2,1);
n_nodes_tot = cell(2,1);
indices = cell(2,1);

totalnodes = 0;
cur = 0;
for i = 1:2
    n_nodes_us{i} = size(fespace_us{i}.nodes,1);
    n_nodes_ps{i} = size(fespace_ps{i}.nodes,1);
    n_nodes_tot{i} = 2*n_nodes_us{i} + n_nodes_ps{i};
    indices{i} = cur+1:cur+2*n_nodes_us{i};
    cur = cur + n_nodes_tot{i};
    totalnodes = totalnodes + n_nodes_tot{i};
end

connectivity = [2 -5;
    0 -4;
    -3 0];

n_boundaries = [5 5];

xc = x1;

B = [];

prods = zeros(3);
for j = 1:3 % interface
    for i = 0:nfreqs(j)
        if (i == 0)
            bx = zeros(totalnodes,1);
            by = zeros(totalnodes,1);
            f = @(x,xp,l) 1;
            for k = 1:2 % domains
                flags = zeros(1,n_boundaries(k));
                si = sign(connectivity(j,k));
                if (si ~= 0)
                    flags(abs(connectivity(j,k))) = si;
                    b1 = zeros(n_nodes_us{k},1);
                    b2 = zeros(n_nodes_us{k},1);
                    b1 = apply_neumann_bc(b1,fespace_us{k},@(x) flags * f(x,xc,ls(j)),8);
                    bx(indices{k},1) = bx(indices{k},1) + [b1;b2];
                    by(indices{k},1) = by(indices{k},1) + [b2;b1];
                end
            end
        else
            bx = zeros(totalnodes,2);
            by = zeros(totalnodes,2);
            f = @(x,xp,l) sin(i*sqrt((x(1,:)-xp(1)).^2+(x(2,:)-xp(2)).^2)/l*pi);
            for k = 1:2 % domains
                flags = zeros(1,n_boundaries(k));
                si = sign(connectivity(j,k));
                if (si ~= 0)
                    flags(abs(connectivity(j,k))) = si;
                    b1 = zeros(n_nodes_us{k},1);
                    b2 = zeros(n_nodes_us{k},1);
                    b1 = apply_neumann_bc(b1,fespace_us{k},@(x) flags * f(x,xc,ls(j)),8);
                    bx(indices{k},1) = bx(indices{k},1) + [b1;b2];
                    by(indices{k},1) = by(indices{k},1) + [b2;b1];
                end
            end
            
            f = @(x,xp,l) cos(i*sqrt((x(1,:)-xp(1)).^2+(x(2,:)-xp(2)).^2)/l*pi);
            for k = 1:2 % domains
                flags = zeros(1,n_boundaries(k));
                si = sign(connectivity(j,k));
                if (si ~= 0)
                    flags(abs(connectivity(j,k))) = si;
                    b1 = zeros(n_nodes_us{k},1);
                    b2 = zeros(n_nodes_us{k},1);
                    b1 = apply_neumann_bc(b1,fespace_us{k},@(x) flags * f(x,xc,ls(j)),8);
                    bx(indices{k},2) = bx(indices{k},2) + [b1;b2];
                    by(indices{k},2) = by(indices{k},2) + [b2;b1];
                end
            end
        end
        B = [B bx by];
    end
end
B = sparse(B);

Bt = B';
dir_indices = [];
cur = 0;
for i = 1:2
    loc_dir_indices = find_dirichlet_indices(fespace_us{i});
    dir_indices = [dir_indices;
        loc_dir_indices + cur;
        loc_dir_indices + cur + n_nodes_us{i}];
    cur = cur + n_nodes_tot{i};
end
B(dir_indices,:) = 0;

nlamb = size(B,2);

jac_block11 = @(u) [];
A = @(u) [];
for i = 1:2
    A = @(u) blkdiag(A(u),As{i}(u(indices{i})));
    jac_block11 = @(u) blkdiag(jac_block11(u),build_jac_navier_stokes(As{i},u(indices{i}),fespace_us{i}));
end

% % solve system with newton's method
% ff = @(u) H(u)*u-rhs;
% tol = 1e-8;
% maxit = 20;
% 
% [sol,er,it] = solve_with_newtons_method(ff,rhs,jac,tol,maxit);
% sols = cell(2,1);
% cur = 0;
% for i = 1:2
%     sols{i}.u1 = sol(cur+1:cur+n_nodes_us{i});
%     sols{i}.u2 = sol(cur+n_nodes_us{i}+1:cur+n_nodes_us{i}*2);
%     sols{i}.p = sol(cur+2*n_nodes_us{i}+1:cur+n_nodes_tot{i});
%     sols{i}.fespace_u = fespace_us{i};
%     sols{i}.fespace_p = fespace_ps{i};
%     cur = cur + n_nodes_tot{i};
% end
% 
% plot_fe_fluid_function(sols{1},'U')
% hold on
% plot_fe_fluid_function(sols{2},'U')

% draw_mesh(mesh_in)
% hold on
% draw_mesh(mesh_out1)
% 
% vertices = [mesh_in.vertices;mesh_out1.vertices;mesh_out2.vertices];
% 
% xm = min(vertices(:,1));
% xM = max(vertices(:,1));
% ym = min(vertices(:,2));
% yM = max(vertices(:,2));
% 
% axis equal
% axis([xm xM ym yM])
% set(gca,'color','none')
% set(gca,'Visible','off')

% ========================= IGA PART =========================

method_data.element_name = 'th';     % Element type for discretization
method_data.degree       = [ 3  3];  % Degree of the splines (pressure space)
method_data.regularity   = [ 2  2];  % Regularity of the splines (pressure space)
method_data.nsub         = [ 6  6];  % Number of subdivisions
method_data.nquad        = [ 5  5];  % Points for the Gaussian quadrature rule

% Penalization parameter for Nitsche's method
factor = 10;
method_data.Cpen = factor*(min(method_data.degree)+1);

[A_iga,B_iga,J_iga,rhs_iga,nurb,int_dofs,vel_iga,Ndofs,space_v,space_p,geometry,msh] = ...
        assemble_system_iga_outflow2(method_data,ls,nfreqs);

B = [B;B_iga];
Bt = [Bt B_iga'];

A = @(u) blkdiag(A(u),A_iga(u(totalnodes+1:end)));
jac_block11 = @(u) blkdiag(jac_block11(u),J_iga(u(totalnodes+1:end)));

H =@(u) [A(u) B;Bt sparse(nlamb,nlamb)];
jac = @(u) [jac_block11(u) B; Bt sparse(nlamb,nlamb)];

x0 = zeros(length(bs{1}) + length(bs{2}) + Ndofs + nlamb,1);

rhs = @(u) [bs{1};bs{2};rhs_iga(u);zeros(nlamb,1)];


% solve system with newton's method
ff = @(u) H(u)*u-rhs(u(totalnodes+1:end));
tol = 1e-8;
maxit = 20;

[sol,er,it] = solve_with_newtons_method(ff,x0,jac,tol,maxit);

sols_fem = cell(2,1);
cur = 0;
for i = 1:2
    sols_fem{i}.u1 = sol(cur+1:cur+n_nodes_us{i});
    sols_fem{i}.u2 = sol(cur+n_nodes_us{i}+1:cur+n_nodes_us{i}*2);
    sols_fem{i}.p = sol(cur+2*n_nodes_us{i}+1:cur+n_nodes_tot{i});
    sols_fem{i}.fespace_u = fespace_us{i};
    sols_fem{i}.fespace_p = fespace_ps{i};
    cur = cur + n_nodes_tot{i};
end

export_vtk_fluid(sols_fem{1},'sol1')
export_vtk_fluid(sols_fem{2},'sol2')

sols_iga = sol(totalnodes+1:end);
v_dofs = sols_iga(1:length(int_dofs));
vel_iga(int_dofs) = v_dofs;

for i = 1:2
    h = trisurf(fespace_us{i}.mesh.elements(:,1:3),fespace_us{i}.mesh.vertices(:,1), ...
        fespace_us{i}.mesh.vertices(:,2),sqrt(sols_fem{i}.u1(1:size(fespace_us{i}.mesh.vertices(:,2),1))).^2 + sols_fem{i}.u2(1:size(fespace_us{i}.mesh.vertices(:,2),1)).^2-10);
    set(h,'edgecolor','none');
    hold on   
end
% shading flat

[eu,ef] = sp_eval(vel_iga,space_v,geometry,[101 101]);
X = squeeze(ef(1,:,:));
Y = squeeze(ef(2,:,:));
Z = squeeze(sqrt(eu(1,:,:).^2 + eu(2,:,:).^2))-10;
surf(X,Y,Z);
view(0,90)
shading interp

axis equal
set(gca,'color','none')
set(gca,'Visible','off')

draw_mesh(mesh_in,[0 0 0])
hold on
draw_mesh(mesh_out1,[0 0 0])

vertices = [mesh_in.vertices;mesh_out1.vertices;mesh_out2.vertices];

vtk_pts = {linspace(0, 1, 51), linspace(0, 1, 51)};

sp_to_vtk (vel_iga, space_v, geometry, vtk_pts, 'sol3', {'velocity', 'divergence'}, {'value', 'divergence'})

xm = min(vertices(:,1));
xM = max(vertices(:,1));
ym = min(vertices(:,2));
yM = max(vertices(:,2));

axis([xm xM ym yM])

hmsh = hierarchical_mesh(msh,[2 2]);
hmsh_plot_cells(hmsh);