function [A,b] = apply_dirichlet_bc_global_matrix_and_rhs(A,b,fespaces,dirichlet_functions)

nfespaces = size(fespaces,2);
dirichlet_indices = [];

aux = 1;
for i = 1:nfespaces
   nnodes = size(fespaces{i}.nodes,1);

   vec = zeros(nnodes,1);
   vec = apply_dirichlet_bc_rhs(vec,fespaces{i},@(x) [1;1;1;1]);
   dirichlet_indices = [dirichlet_indices;find(vec == 1)+aux-1];
   
   b(aux:aux+nnodes-1) = apply_dirichlet_bc_rhs(b(aux:aux+nnodes-1), ...
                            fespaces{i},dirichlet_functions);
   aux = aux + nnodes;
end
n_indices = length(dirichlet_indices);

A(dirichlet_indices,:) = 0;
A(dirichlet_indices,dirichlet_indices) = spdiags(ones(n_indices,1)*1, ...
                                            0,n_indices,n_indices);