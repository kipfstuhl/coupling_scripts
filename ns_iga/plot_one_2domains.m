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

ref_inflow = 11;
ref_out1 = 4;

% ========================= FINITE ELEMENT PART =========================

% load finite element matrices (we load the second outflow just to retrieve
% geometric information about the interface 3)
mesh_in = read_mesh(['../meshes_old/refinement',num2str(ref_inflow),'/inflow_distorted.msh']);
mesh_out1 = read_mesh(['../meshes_old/refinement',num2str(ref_out1),'/outflow1_distorted.msh']);
mesh_out2 = read_mesh(['../meshes_old/refinement',num2str(ref_out1),'/outflow2_distorted.msh']);
h = max(mesh_in.h,mesh_out1.h);

fespace_us = cell(2,1);
fespace_ps = cell(2,1);

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

nfreqs = [5 5 5];

n_nodes_us = cell(2,1);
n_nodes_ps = cell(2,1);
n_nodes_tot = cell(2,1);
indices = cell(2,1);

totalnodes = 0;
cur = 0;
for i = 2:2
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
            for k = 2:2 % domains
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
            for k = 2:2 % domains
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
            for k = 2:2 % domains
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
save('B.mat','B')

Bt = B';
dir_indices = [];
cur = 0;
for i = 2:2
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
for i = 2:2
    A = @(u) blkdiag(A(u),As{i}(u(indices{i})));
    jac_block11 = @(u) blkdiag(jac_block11(u),build_jac_navier_stokes(As{i},u(indices{i}),fespace_us{i}));
end

A_iga = cell(2,1);
B_iga = cell(2,1);
J_iga = cell(2,1);
rhs_iga = cell(2,1);
nurb = cell(2,1);
int_dofs = cell(2,1);
vel_iga = cell(2,1);
Ndofs = cell(2,1);
space_v = cell(2,1);
space_p = cell(2,1);
geometry = cell(2,1);
msh = cell(2,1);
dofs_iga = cell(2,1);
% ========================= IGA PART INFLOW =========================

method_data.element_name = 'th';     % Element type for discretization
method_data.degree       = [ 3  3];  % Degree of the splines (pressure space)
method_data.regularity   = [ 2  2];  % Regularity of the splines (pressure space)
method_data.nsub         = [ 6  6];  % Number of subdivisions
method_data.nquad        = [ 5  5];  % Points for the Gaussian quadrature rule

% Penalization parameter for Nitsche's method
factor = 10;
method_data.Cpen = factor*(min(method_data.degree)+1);

ii = 1;
[A_iga{ii},B_iga{ii},J_iga{ii},rhs_iga{ii},nurb{ii},int_dofs{ii},vel_iga{ii},Ndofs{ii},space_v{ii},space_p{ii},geometry{ii},msh{ii}] = ...
    assemble_system_iga_inflow(method_data,ls,nfreqs);

% ========================= IGA PART OUTFLOW 2 =========================

method_data.element_name = 'th';     % Element type for discretization
method_data.degree       = [ 2  2];  % Degree of the splines (pressure space)
method_data.regularity   = [ 1  1];  % Regularity of the splines (pressure space)
method_data.nsub         = [ 5  5];  % Number of subdivisions
method_data.nquad        = [ 5  5];  % Points for the Gaussian quadrature rule

% Penalization parameter for Nitsche's method
factor = 10;
method_data.Cpen = factor*(min(method_data.degree)+1);

ii = 2;
[A_iga{ii},B_iga{ii},J_iga{ii},rhs_iga{ii},nurb{ii},int_dofs{ii},vel_iga{ii},Ndofs{ii},space_v{ii},space_p{ii},geometry{ii},msh{ii}] = ...
    assemble_system_iga_outflow2(method_data,ls,nfreqs);

aux = totalnodes;
for i = 1:2
    B = [B;B_iga{i}];
    Bt = [Bt B_iga{i}'];
    dofs_iga{i} = aux+1:aux+Ndofs{i};
    A = @(u) blkdiag(A(u),A_iga{i}(u(dofs_iga{i})));
    jac_block11 = @(u) blkdiag(jac_block11(u),J_iga{i}(u(dofs_iga{i})));
    aux = aux+Ndofs{i};
end

H =@(u) [A(u) B;Bt sparse(nlamb,nlamb)];
jac = @(u) [jac_block11(u) B; Bt sparse(nlamb,nlamb)];

x0 = zeros(length(bs{2}) + Ndofs{1} + Ndofs{2} + nlamb,1);

rhs = @(u) [bs{2};rhs_iga{1}(u(dofs_iga{1}));rhs_iga{2}(u(dofs_iga{2}));zeros(nlamb,1)];

% solve system with newton's method
ff = @(u) H(u)*u-rhs(u);
tol = 1e-8;
maxit = 20;

[sol,er,it] = solve_with_newtons_method(ff,x0,jac,tol,maxit);

sols_fem = cell(2,1);
cur = 0;
for i = 2:2
    sols_fem{i}.u1 = sol(cur+1:cur+n_nodes_us{i});
    sols_fem{i}.u2 = sol(cur+n_nodes_us{i}+1:cur+n_nodes_us{i}*2);
    sols_fem{i}.p = sol(cur+2*n_nodes_us{i}+1:cur+n_nodes_tot{i});
    sols_fem{i}.fespace_u = fespace_us{i};
    sols_fem{i}.fespace_p = fespace_ps{i};
    cur = cur + n_nodes_tot{i};
end

sol_iga = sol(dofs_iga{1});
v_dofs = sol_iga(1:length(int_dofs{1}));
vel_iga{1}(int_dofs{1}) = v_dofs;
p_iga{1} = sol_iga(length(int_dofs{1})+1:end);


sol_iga = sol(dofs_iga{2});
v_dofs = sol_iga(1:length(int_dofs{2}));
vel_iga{2}(int_dofs{2}) = v_dofs;
p_iga{2} = sol_iga(length(int_dofs{2})+1:end);


% ========================= COMPUTING ERROR =========================

%         Au = op_gradu_gradv_tp (space_v, space_v, msh, @(x, y) ones (size (x)));
%         Mu = op_u_v_tp (space_v, space_v, msh, @(x, y) ones (size (x)));
%
%         energy(1) = sqrt(vel_iga'*(Au+Mu)*vel_iga);
%
%         for i = 1:2
%             Ai = assemble_stiffness(1,fespace_us{i});
%             Mi = assemble_mass(fespace_us{i});
%             matnorm = blkdiag(Ai+Mi,Ai+Mi);
%             uss = [sols_fem{i}.u1;sols_fem{i}.u2];
%             energy(i+1) = sqrt(uss'*matnorm*uss);
%         end
%
%         err_v(count,countdeg) = sqrt(abs(energy*energy'-exnorm_u^2))

i = 2;
h = trisurf(fespace_us{i}.mesh.elements(:,1:3),fespace_us{i}.mesh.vertices(:,1), ...
    fespace_us{i}.mesh.vertices(:,2),sqrt(sols_fem{i}.u1(1:size(fespace_us{i}.mesh.vertices(:,2),1))).^2 + sols_fem{i}.u2(1:size(fespace_us{i}.mesh.vertices(:,2),1)).^2-10);
set(h,'edgecolor','none');
hold on
% shading flat

for i = 1:2
    [eu,ef] = sp_eval(vel_iga{i},space_v{i},geometry{i},[101 101]);
    X = squeeze(ef(1,:,:));
    Y = squeeze(ef(2,:,:));
    Z = squeeze(sqrt(eu(1,:,:).^2 + eu(2,:,:).^2))-10;
    surf(X,Y,Z);
    view(0,90)
end
shading interp

axis equal
set(gca,'color','none')
set(gca,'Visible','off')

draw_mesh(mesh_out1,[0 0 0])
hold on

vertices = [mesh_in.vertices;mesh_out1.vertices;mesh_out2.vertices];

vtk_pts = {linspace(0, 1, 51), linspace(0, 1, 51)};

% sp_to_vtk (vel_iga{1}, space_v{2}, geometry{2}, vtk_pts, 'sol3', {'velocity', 'divergence'}, {'value', 'divergence'})
% 
% sp_to_vtk (vel_iga{2}, space_v{2}, geometry{2}, vtk_pts, 'sol2', {'velocity', 'divergence'}, {'value', 'divergence'})

for i = 1:2
    hmsh = hierarchical_mesh(msh{i},[2 2]);
    hmsh_plot_cells(hmsh);   
end

xm = min(vertices(:,1));
xM = max(vertices(:,1));
ym = min(vertices(:,2));
yM = max(vertices(:,2));

axis([xm xM ym yM])

hold on

create_inflow_nurb
for i = 1:ni
    for j = 1:nj
        p = coefs(:,i,j);
        plot(p(1),p(2),'.r','Markersize',20);
    end
end

create_outflow2_nurb
for i = 1:ni
    for j = 1:nj
        p = coefs(:,i,j);
        plot(p(1),p(2),'.r','Markersize',20);
    end
end
