function [blocks,blocks_t] = couple_navier_stokes_solutions_single_basis(fespaces_u,fespaces_p,bcs_flags,bf,gausspoints)

ncouplings = size(fespaces_u,2);

blocks_x = [];
blocks_y = [];
blocks_xt = [];
blocks_yt = [];
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
    bsx{j} = apply_neumann_bc(bsx{j},fespaces_u{j},@(x) bf(x)*bcs_flags(:,j),gausspoints);
    blocks_x = [blocks_x [bsx{j};bsy{j};bsp{j}]'];
    blocks_y = [blocks_y [bsy{j};bsx{j};bsp{j}]'];
    
    bsx{j} = apply_dirichlet_bc_rhs(bsx{j},fespaces_u{j},@(x) [0;0;0;0]);
    blocks_xt = [blocks_xt;[bsx{j};bsy{j};bsp{j}]];
    blocks_yt = [blocks_yt;[bsy{j};bsx{j};bsp{j}]];
end

blocks = sparse([blocks_x;blocks_y]);
blocks_t = sparse([blocks_xt blocks_yt]);

