clear all
close all
clc

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex')

ref_inflow = 1;
ref_out1 = 1;
ref_out2 = 1;

mesh_in = read_mesh(['../meshes/refinement',num2str(ref_inflow),'/inflow.msh']);
mesh_out1 = read_mesh(['../meshes/refinement',num2str(ref_out1),'/outflow1.msh']);
mesh_out2 = read_mesh(['../meshes/refinement',num2str(ref_out2),'/outflow2.msh']);

fespace_us = cell(3,1);
fespace_ps = cell(3,1);

% fespaces inlet
fespace_us{1} = create_fespace(mesh_in,'P2',[1 0 0 1 1]);
fespace_ps{1} = create_fespace(mesh_in,'P1',[]);

% fespaces outlet1
fespace_us{2} = create_fespace(mesh_out1,'P2',[1 0 1 0 0]);
fespace_ps{2} = create_fespace(mesh_out1,'P1',[]);

% fespaces outlet2
fespace_us{3} = create_fespace(mesh_out2,'P2',[0 1 0 1 0]);
fespace_ps{3} = create_fespace(mesh_out2,'P1',[]);

r_in = 0.5;
U = 10;
dir0 = @(x) zeros(2,6);
dir_in = @(x) [0 0; 0 0; 0 0;0 0; (r_in^2 - x(2)^2)/r_in^2 * U 0]';
neu0 = dir0;

As = cell(3,1);
bs = cell(3,1);

[As{1},bs{1}] = assembler_steady_navier_stokes(fespace_us{1},fespace_ps{1},[0;0],1,dir_in,neu0);
[As{2},bs{2}] = assembler_steady_navier_stokes(fespace_us{2},fespace_ps{2},[0;0],1,dir0,neu0);
[As{3},bs{3}] = assembler_steady_navier_stokes(fespace_us{3},fespace_ps{3},[0;0],1,dir0,neu0);

nfreq = 10;
nfourier = nfreq*2+1;

% find interface 1
b1 = mesh_in.boundaries{2};
x1 = mesh_in.vertices(b1(1,1),1:2);
x2 = mesh_in.vertices(b1(end,end),1:2);
dif = x2 - x1;
l1 = norm(dif);
t1 = dif/l1;
n1 = [t1(2) -t1(1)];

x1
x2

% find interface 2
b2 = mesh_out1.boundaries{4};
x2 = mesh_out1.vertices(b2(1,1),1:2);
x1 = mesh_out1.vertices(b2(end,end),1:2);
dif = x2 - x1;
l2 = norm(dif);
t2= dif/l2;
n2 = [t2(2) -t2(1)];

plot(x1(1),x1(2),'.','Markersize',20)
plot(x2(1),x2(2),'.','Markersize',20)


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

% build mass matrices on sample points
% nsamples = 1000;
% 
% R1 = compute_orthonormalization_matrix_onlysin(nfreq,nsamples,ls(1),0.5);
% R2 = compute_orthonormalization_matrix_onlysin(nfreq,nsamples,ls(2),0.5);
% R3 = compute_orthonormalization_matrix_onlysin(nfreq,nsamples,ls(3),0.5);
% R = blkdiag(R1,R1,R2,R2,R3,R3);

n_nodes_us = cell(3,1);
n_nodes_ps = cell(3,1);
n_nodes_tot = cell(3,1);
indices = cell(3,1);
totalnodes = 0;
cur = 0;
for i = 1:3
    n_nodes_us{i} = size(fespace_us{i}.nodes,1);
    n_nodes_ps{i} = size(fespace_ps{i}.nodes,1);
    n_nodes_tot{i} = 2*n_nodes_us{i} + n_nodes_ps{i};
    indices{i} = cur+1:cur+2*n_nodes_us{i};
    cur = cur + n_nodes_tot{i};
    totalnodes = totalnodes + n_nodes_tot{i};
end

connectivity = [2 -5 0;
    0 4 -1;
    3 0 -5];

n_boundaries = [5 5 5];

xc = x1;

B = [];

% add "constants"
for i = 1:3 % interface 1
    bx = zeros(totalnodes,1);
    by = zeros(totalnodes,1);
    for j = 1:3 % interface 2
        prod = abs(ts(i,:)*ns(j,:)');
        disp(['i = ', num2str(i),' j = ',num2str(j),' prod = ',num2str(prod)])
        for k = 1:3 % domains
            flags = zeros(1,n_boundaries(k));
            si = sign(connectivity(j,k));
            if (si ~= 0)
                flags(abs(connectivity(j,k))) = si;
                b1 = zeros(n_nodes_us{k},1);
                b2 = zeros(n_nodes_us{k},1);
                b1 = apply_neumann_bc(b1,fespace_us{k},@(x) prod*flags,8);
                bx(indices{k}) = bx(indices{k}) + [b1;b2];
                by(indices{k}) = by(indices{k}) + [b2;b1];
            end
        end
    end
    B = [B bx by];
end

for j = 1:3 % interface
    for i = 1:nfreq
        bx = zeros(totalnodes,1);
        by = zeros(totalnodes,1);
        if (i == 0)
            f = @(x,xp,l) 1;
        else
            f = @(x,xp,l) sin(i*sqrt((x(1,:)-xp(1)).^2+(x(2,:)-xp(2)).^2)/l*pi/2);
        end
        for k = 1:3 % domains
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
        
        %         f = @(x,xp,l) cos(i*norm(x-xp)/l*pi);
        %         for k = 1:3 % domains
        %             flags = zeros(1,n_boundaries(k));
        %             si = sign(connectivity(j,k));
        %             if (si ~= 0)
        %                 flags(abs(connectivity(j,k))) = si;
        %                 b1 = zeros(n_nodes_us{k},1);
        %                 b2 = zeros(n_nodes_us{k},1);
        %                 b1 = apply_neumann_bc(b1,fespace_us{k},@(x) flags * f(x,xc,ls(j)),8);
        %                 bx(indices{k},2) = bx(indices{k},2) + [b1;b2];
        %                 by(indices{k},2) = by(indices{k},2) + [b2;b1];
        %             end
        %         end
        B = [B bx by];
    end
end
B = sparse(B);

Bt = B';
% Bt = R'\Bt;
B = Bt';
dir_indices = [];
cur = 0;
for i = 1:3
    loc_dir_indices = find_dirichlet_indices(fespace_us{i});
    dir_indices = [dir_indices;
        loc_dir_indices + cur;
        loc_dir_indices + cur + n_nodes_us{i}];
%     for j = loc_dir_indices(:)'
%         plot(fespace_us{i}.nodes(j,1),fespace_us{i}.nodes(j,2),'.r','Markersize',10);
%     end
    cur = cur + n_nodes_tot{i};
end
B(dir_indices,:) = 0;

nlamb = size(B,2);

A = @(u) [As{1}(u) sparse(n_nodes_tot{1},n_nodes_tot{2}) sparse(n_nodes_tot{1},n_nodes_tot{3}); ...
    sparse(n_nodes_tot{2},n_nodes_tot{1}) As{2}(u) sparse(n_nodes_tot{2},n_nodes_tot{3}); ...
    sparse(n_nodes_tot{3},n_nodes_tot{1}) sparse(n_nodes_tot{3},n_nodes_tot{2}) As{3}(u)];

jac_block11 = @(u) [];
for i = 1:3
    jac_block11 = @(u) blkdiag(jac_block11(u),build_jac_navier_stokes(As{i},u(indices{i}),fespace_us{i}));
end

jac = @(u) [jac_block11(u) B; Bt sparse(nlamb,nlamb)];

H =@(u) [A(u) B;Bt sparse(nlamb,nlamb)];
rhs = [bs{1};bs{2};bs{3};zeros(nlamb,1)];

u0 = zeros(length(rhs),1);
x0 = H(u0)\rhs;

% solve system with newton's method
f = @(u) H(u)*u-rhs;
tol = 1e-8;
maxit = 20;

[sol,er,it] = solve_with_newtons_method(f,x0,jac,tol,maxit);

sols = cell(3,1);
cur = 0;
for i = 1:3
    sols{i}.u1 = sol(cur+1:cur+n_nodes_us{i});
    sols{i}.u2 = sol(cur+n_nodes_us{i}+1:cur+n_nodes_us{i}*2);
    sols{i}.p = sol(cur+2*n_nodes_us{i}+1:cur+n_nodes_tot{i});
    sols{i}.fespace_u = fespace_us{i};
    sols{i}.fespace_p = fespace_ps{i};
    cur = cur + n_nodes_tot{i};
end

plot_fe_fluid_function(sols{1},'U')
hold on
plot_fe_fluid_function(sols{2},'U')
plot_fe_fluid_function(sols{3},'U')
draw_mesh(mesh_in,[0 0 0]);
draw_mesh(mesh_out1,[0 0 0]);
draw_mesh(mesh_out2,[0 0 0]);

errsu = zeros(3,1);
errsp = zeros(3,1);

load('exsol.mat');
exsol = sol;

exsol_u = @(x) [evaluate_fe_function(exsol.u1,exsol.fespace_u,x); ...
                evaluate_fe_function(exsol.u2,exsol.fespace_u,x)];
exsol_p = @(x) evaluate_fe_function(exsol.p,exsol.fespace_p,x);

exsol_grad = @(x) [evaluate_fe_function_gradient(exsol.u1,exsol.fespace_u,x)'; ...
                   evaluate_fe_function_gradient(exsol.u2,exsol.fespace_u,x)'];
for i = 1:3
    errsu(i) = compute_H1_error_velocity(fespace_us{i},sols{i},exsol_u,exsol_grad);
    errsp(i) = compute_L2_error(fespace_ps{i},sols{i}.p,exsol_p);
end

errs = errsu + errsp;
err = sqrt(errs'*errs);

keyboard

vertices = [mesh_in.vertices;mesh_out1.vertices;mesh_out2.vertices];

xm = min(vertices(:,1));
xM = max(vertices(:,1));
ym = min(vertices(:,2));
yM = max(vertices(:,2));

axis([xm xM ym yM])

% plot_fe_function(sols{1}.u1,fespace_us{1})
% plot_fe_function(sols{2}.u1,fespace_us{2})
% plot_fe_function(sols{3}.u1,fespace_us{3})

set(gca,'color','none') 
set(gca,'Visible','off') 
export_vtk_fluid(sols{1},'sol1')
export_vtk_fluid(sols{2},'sol2')
export_vtk_fluid(sols{3},'sol3')
