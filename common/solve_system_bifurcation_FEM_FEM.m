function [sols,h] = solve_system_bifurcation_FEM_FEM(ref,freqq)

meshes = cell(3,1);
mesh_labels = {'inflow_distorted.msh', ...
    'outflow1_distorted.msh', ...
    'outflow2_distorted.msh'};
h = 0;

% read mesh
for i = 1:3
    meshes{i} = read_mesh(['../meshes/refinement',num2str(ref),'/',mesh_labels{i}]);
    h = max(h,meshes{i}.h);
end

fespace_us = cell(3,1);
fespace_ps = cell(3,1);

% fespaces inlet
fespace_us{1} = create_fespace(meshes{1},'P2',[1 0 0 1 1]);
fespace_ps{1} = create_fespace(meshes{1},'P1',[]);

% fespaces outlet1
fespace_us{2} = create_fespace(meshes{2},'P2',[1 0 1 0 0]);
fespace_ps{2} = create_fespace(meshes{2},'P1',[]);

% fespaces outlet2
fespace_us{3} = create_fespace(meshes{3},'P2',[0 1 0 1 0]);
fespace_ps{3} = create_fespace(meshes{3},'P1',[]);

% parameters for the boundary conditions
r_in = 0.5;
U = 1;
dir0 = @(x) zeros(2,6);
dir_in = @(x) [0 0; 0 0; 0 0;0 0; (r_in^2 - x(2)^2)/r_in^2 * U 0]';
neu0 = dir0;

% build matrices and rhs
As = cell(3,1);
bs = cell(3,1);

[As{1},bs{1}] = assembler_steady_navier_stokes(fespace_us{1},fespace_ps{1},[0;0],1,dir_in,neu0);
[As{2},bs{2}] = assembler_steady_navier_stokes(fespace_us{2},fespace_ps{2},[0;0],1,dir0,neu0);
[As{3},bs{3}] = assembler_steady_navier_stokes(fespace_us{3},fespace_ps{3},[0;0],1,dir0,neu0);

% find interface 1 and compute length
b1 = meshes{1}.boundaries{2};
x1 = meshes{1}.vertices(b1(end,end),1:2);
x2 = meshes{1}.vertices(b1(1,1),1:2);
dif = x2 - x1;
l1 = norm(dif);

% find interface 2 and compute length
b2 = meshes{2}.boundaries{4};
x1 = meshes{2}.vertices(b2(1,1),1:2);
x2 = meshes{2}.vertices(b2(end,end),1:2);
dif = x2 - x1;
l2 = norm(dif);

% find interface 3 and compute length
b3 = meshes{3}.boundaries{5};
x1 = meshes{3}.vertices(b3(1,1),1:2);
x2 = meshes{3}.vertices(b3(end,end),1:2);
dif = x2 - x1;
l3 = norm(dif);

ls = [l1 l2 l3];

nfreqs = [1 1 1] * freqq;

% store dimensions of the spaces and offsets
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

% connectivity structure of the domain. The component (i,j) contains the 0
% if the jth domain does not communicate with interface i, otherwise it
% contains the flag of the corresponding interface (positive sign = inward
% normal, negative sign = outward normal)
connectivity = [2 -5 0;
    0 4 -1;
    -3 0 5];

% number of boundaries (i.e. flags) for each subdomain
n_boundaries = [5 5 5];

% cross point
xc = x1;

B = [];

prods = zeros(3);

% we manually build the coupling matrices
for j = 1:3 % interface
    for i = 0:nfreqs(j) % frequencies
        if (i == 0) % => then we integrate just the constant
            bx = zeros(totalnodes,1);
            by = zeros(totalnodes,1);
            f = @(x,xp,l) 1;
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
        else
            bx = zeros(totalnodes,2);
            by = zeros(totalnodes,2);
            
            f = @(x,xp,l) sin(i*sqrt((x(1,:)-xp(1)).^2+(x(2,:)-xp(2)).^2)/l*pi);
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
            
            f = @(x,xp,l) cos(i*sqrt((x(1,:)-xp(1)).^2+(x(2,:)-xp(2)).^2)/l*pi);
            for k = 1:3 % domains
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

% save Bt
Bt = B';

% and now apply bcs to B
dir_indices = [];
cur = 0;
for i = 1:3
    loc_dir_indices = find_dirichlet_indices(fespace_us{i});
    dir_indices = [dir_indices;
        loc_dir_indices + cur;
        loc_dir_indices + cur + n_nodes_us{i}];
    cur = cur + n_nodes_tot{i};
end
B(dir_indices,:) = 0;

nlamb = size(B,2);

% build the global system + jacobian
jac_block11 = @(u) [];
A = @(u) [];
for i = 1:3
    A = @(u) blkdiag(A(u),As{i}(u(indices{i})));
    jac_block11 = @(u) blkdiag(jac_block11(u),build_jac_navier_stokes(As{i},u(indices{i}),fespace_us{i}));
end


% solve the non linear system
x0 = zeros(length(bs{1}) + length(bs{2}) + length(bs{3}),1);
x0 = [x0;zeros(nlamb,1)];
jac = @(u) [jac_block11(u) B; Bt sparse(nlamb,nlamb)];

H =@(u) [A(u) B;Bt sparse(nlamb,nlamb)];

rhs = [bs{1};bs{2};bs{3};zeros(nlamb,1)];

% solve system with newton's method
ff = @(u) H(u)*u-rhs;
tol = 1e-8;
maxit = 20;

[sol,~,~] = solve_with_newtons_method(ff,x0,jac,tol,maxit);

% separe the solutions into among the subdomains
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