function [M,R] = build_interface_matrix_and_restriction(fespace,boundary_index)
if (~strcmp(fespace.degree,'P2'))
    error(['Interface matrix not implemented for ',fespace.degree,' elements!']);
end
mesh = fespace.mesh;
nodes = fespace.nodes;

nnodes = size(nodes,1);
indicesonboundary = [];
for i = 1:nnodes
    if (nodes(i,3) == boundary_index || nodes(i,4) == boundary_index)
        indicesonboundary = [indicesonboundary;i];
    end
end

nindices = length(indicesonboundary);

if (boundary_index == 2 || boundary_index == 4)
    M = build_1D_mass(mesh.yp,mesh.yp+mesh.H,mesh.h2);
else
    M = build_1D_mass(mesh.xp,mesh.xp+mesh.L,mesh.h1);
end

R = zeros(nindices,nnodes);
for i = 1:nindices
    R(i,indicesonboundary(i)) = 1;
end
R = sparse(R);