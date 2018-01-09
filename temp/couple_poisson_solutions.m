function [blocks,blocks_t] = couple_poisson_solutions(fespaces,bcs_flags,bf,gausspoints)

ncouplings = size(fespaces,2);

blocks = [];
blocks_t = [];
nnodes = zeros(ncouplings,1);

for i = 1:ncouplings
    nnodes(i) = size(fespaces{i}.nodes,1);
end

bsx = {};

for j = 1:ncouplings
    bsx{end+1} = zeros(nnodes(j),1);
end

for j = 1:ncouplings
    bsx{j} = apply_neumann_bc(bsx{j},fespaces{j},@(x) bf(x)*bcs_flags(:,j),gausspoints);
    blocks = [blocks bsx{j}'];
    
    bsx{j} = apply_dirichlet_bc_rhs(bsx{j},fespaces{j},@(x) [0;0;0;0]);
    blocks_t = [blocks_t;bsx{j}];
end

blocks = sparse(blocks);
blocks_t = sparse(blocks_t);



