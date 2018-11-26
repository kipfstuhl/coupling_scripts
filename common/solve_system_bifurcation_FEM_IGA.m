function [sol_fem,is] = solve_system_bifurcation_FEM_IGA(options_inflow,options_outflow1,options_outflow2,freq)

% we start by assembling the structures for the finite element method. For
% the sake of clarity, we do not stick to the order of the domain, but we 
% assemble first the fem domain and then the iga domains
mesh_in = read_mesh(['../meshes/refinement',num2str(options_outflow1.ref),'_fem_iga/inflow_distorted.msh']);
mesh_out1 = read_mesh(['../meshes/refinement',num2str(options_outflow1.ref),'_fem_iga/outflow1_distorted.msh']);
mesh_out2 = read_mesh(['../meshes/refinement',num2str(options_outflow1.ref),'_fem_iga/outflow2_distorted.msh']);

% fespaces outlet1
fespace_u = create_fespace(mesh_out1,options_outflow1.polydegree_u,[1 0 1 0 0]);
fespace_p = create_fespace(mesh_out1,options_outflow1.polydegree_p,[]);

% boundary conditions
dir0 = @(x) zeros(2,6);
neu0 = dir0;

[As,bs] = assembler_steady_navier_stokes(fespace_u,fespace_p,[0;0],1,dir0,neu0);

% find interface 1 and compute length
b1 = mesh_in.boundaries{2};
x1 = mesh_in.vertices(b1(end,end),1:2);
x2 = mesh_in.vertices(b1(1,1),1:2);
dif = x2 - x1;
l1 = norm(dif);

% find interface 2 and compute length
b2 = mesh_out1.boundaries{4};
x1 = mesh_out1.vertices(b2(1,1),1:2);
x2 = mesh_out1.vertices(b2(end,end),1:2);
dif = x2 - x1;
l2 = norm(dif);

% find interface 3 and compute length
b3 = mesh_out2.boundaries{5};
x1 = mesh_out2.vertices(b3(1,1),1:2);
x2 = mesh_out2.vertices(b3(end,end),1:2);
dif = x2 - x1;
l3 = norm(dif);

ls = [l1 l2 l3];

nfreqs = [1 1 1]*freq;

connectivity = [-5;-4];

n_boundaries = 5;

xc = x1;

B = [];

n_nodes_u = size(fespace_u.nodes,1);
n_nodes_p = size(fespace_p.nodes,1);
totalnodes = 2*n_nodes_u + n_nodes_p;

% manually assemble coupling part for finite element space
for j = 1:2 % interface
    for i = 0:nfreqs(j)
        if (i == 0)
            bx = zeros(totalnodes,1);
            by = zeros(totalnodes,1);
            f = @(x,xp,l) 1;
            flags = zeros(1,n_boundaries);
            si = sign(connectivity(j));
            if (si ~= 0)
                flags(abs(connectivity(j))) = si;
                b1 = zeros(n_nodes_u,1);
                b2 = zeros(n_nodes_u,1);
                b1 = apply_neumann_bc(b1,fespace_u,@(x) flags * f(x,xc,ls(j)),8);
                bx(1:n_nodes_u*2) = bx(1:n_nodes_u*2) + [b1;b2];
                by(1:n_nodes_u*2) = by(1:n_nodes_u*2) + [b2;b1];
            end
        else
            bx = zeros(totalnodes,2);
            by = zeros(totalnodes,2);
            f = @(x,xp,l) sin(i*sqrt((x(1,:)-xp(1)).^2+(x(2,:)-xp(2)).^2)/l*pi);
            flags = zeros(1,n_boundaries);
            si = sign(connectivity(j));
            if (si ~= 0)
                flags(abs(connectivity(j))) = si;
                b1 = zeros(n_nodes_u,1);
                b2 = zeros(n_nodes_u,1);
                b1 = apply_neumann_bc(b1,fespace_u,@(x) flags * f(x,xc,ls(j)),8);
                bx(1:n_nodes_u*2,1) = bx(1:n_nodes_u*2,1) + [b1;b2];
                by(1:n_nodes_u*2,1) = by(1:n_nodes_u*2,1) + [b2;b1];
            end
            
            f = @(x,xp,l) cos(i*sqrt((x(1,:)-xp(1)).^2+(x(2,:)-xp(2)).^2)/l*pi);
            flags = zeros(1,n_boundaries);
            si = sign(connectivity(j));
            if (si ~= 0)
                flags(abs(connectivity(j))) = si;
                b1 = zeros(n_nodes_u,1);
                b2 = zeros(n_nodes_u,1);
                b1 = apply_neumann_bc(b1,fespace_u,@(x) flags * f(x,xc,ls(j)),8);
                bx(1:n_nodes_u*2,2) = bx(1:n_nodes_u*2,2) + [b1;b2];
                by(1:n_nodes_u*2,2) = by(1:n_nodes_u*2,2) + [b2;b1];
            end
        end
        B = [B bx by];
    end
end
% we add zeros for the part corresponding to the missing interface
Bzeros = zeros(size(B,1),size(B,2)/2);
B = sparse([B Bzeros]);
Bt = B';

% apply boundary conditions
dir_indices = find_dirichlet_indices(fespace_u);
B([dir_indices;dir_indices+n_nodes_u],:) = 0;

nlamb = size(B,2);

jac_block11 = @(u) [];

% start assembling matrices
A = @(u) As(u(1:n_nodes_u*2));
jac_block11 = @(u) blkdiag(jac_block11(u),build_jac_navier_stokes(As,u(1:n_nodes_u*2),fespace_u));

% iga part
is{1} = assemble_system_iga_inflow(options_inflow,ls,freq*[1 1 1]);
is{2} = assemble_system_iga_outflow2(options_outflow2,ls,freq*[1 1 1]);

aux = totalnodes;
for i = 1:2
    B = [B; is{i}.B];
    Bt = [Bt is{i}.B'];
    is{i}.dofs_iga = aux+1:aux+is{i}.Ndofs;
    A = @(u) blkdiag(A(u),is{i}.A(u(is{i}.dofs_iga)));
    jac_block11 = @(u) blkdiag(jac_block11(u),is{i}.J(u(is{i}.dofs_iga)));
    aux = aux+is{i}.Ndofs;
end

H =@(u) [A(u) B;Bt sparse(nlamb,nlamb)];
jac = @(u) [jac_block11(u) B; Bt sparse(nlamb,nlamb)];

x0 = zeros(length(bs) + is{1}.Ndofs + is{2}.Ndofs + nlamb,1);

rhs = @(u) [bs; ...
            is{1}.rhs(u(is{1}.dofs_iga)); ...
            is{2}.rhs(u(is{2}.dofs_iga)); ...
            zeros(nlamb,1)];

% solve system with newton's method
ff = @(u) H(u)*u-rhs(u);
tol = 1e-8;
maxit = 20;

[sol,~,~] = solve_with_newtons_method(ff,x0,jac,tol,maxit);

cur = 0;

% retrieve fem solution
sol_fem.u1 = sol(cur+1:cur+n_nodes_u);
sol_fem.u2 = sol(cur+n_nodes_u+1:cur+n_nodes_u*2);
sol_fem.p = sol(cur+2*n_nodes_u+1:cur+totalnodes);
sol_fem.fespace_u = fespace_u;
sol_fem.fespace_p = fespace_p;

% retrieve iga solutions
for i = 1:2
    sol_iga = sol(is{i}.dofs_iga);
    v_dofs = sol_iga(1:length(is{i}.int_dofs));
    is{i}.vel(is{i}.int_dofs) = v_dofs;
    p = sol_iga(length(is{i}.int_dofs)+1:end);
    is{i}.p = p;
end
