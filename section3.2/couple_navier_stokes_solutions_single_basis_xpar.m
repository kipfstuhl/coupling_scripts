function [blocks,blocks_t] = couple_navier_stokes_solutions_single_basis_xpar(fespaces_u,fespaces_p,bcs_flags,bf,gausspoints)

ncouplings = size(fespaces_u,2);

blocks_xx = [];
blocks_xy = [];
blocks_yx = [];
blocks_yy = [];

blocks_xxt = [];
blocks_xyt = [];
blocks_yxt = [];
blocks_yyt = [];

nnodes_u = zeros(ncouplings,1);
nnodes_p = zeros(ncouplings,1);

for i = 1:ncouplings
    nnodes_u(i) = size(fespaces_u{i}.nodes,1);
    nnodes_p(i) = size(fespaces_p{i}.nodes,1);
end

bsx = {};
bsy = {};
bsp = {};

for j = 1:ncouplings
    bsx{end+1} = zeros(nnodes_u(j),1);
    bsy{end+1} = zeros(nnodes_u(j),1);
    bsp{end+1} = zeros(nnodes_p(j),1);
end

for j = 1:ncouplings    
    bc_flags_aux = [bcs_flags(1,j) 0 bcs_flags(3,j) 0]';
    bsx{j} = apply_neumann_bc(bsx{j},fespaces_u{j},@(x) bf(x)*bc_flags_aux,gausspoints);
    blocks_xy = [blocks_xy [bsx{j};bsy{j};bsp{j}]'];
    blocks_yy = [blocks_yy [bsy{j};bsx{j};bsp{j}]'];
    
    bsx{j} = apply_dirichlet_bc_rhs(bsx{j},fespaces_u{j},@(x) [0;0;0;0]);
    blocks_xyt = [blocks_xyt;[bsx{j};bsy{j};bsp{j}]];
    blocks_yyt = [blocks_yyt;[bsy{j};bsx{j};bsp{j}]];
end

blocks = sparse([blocks_xx;blocks_xy;blocks_yx;blocks_yy]);
blocks_t = sparse([blocks_xxt blocks_xyt blocks_yxt blocks_yyt]);

