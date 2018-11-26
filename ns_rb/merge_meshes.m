function mesh = merge_meshes(mesh1,mesh2,bindex1,bindex2,boundary_map)

bindices1 = [mesh1.boundaries{bindex1}(:,1);mesh1.boundaries{bindex1}(:,end)];
bindices2 = [mesh2.boundaries{bindex2}(:,1);mesh2.boundaries{bindex2}(:,end)];

if (length(bindices1) ~= length(bindices2))
    error('Meshes are not mergeable!');
end

if (norm(mesh1.vertices(bindices1(1),1:2) - mesh2.vertices(bindices2(1),1:2)) > 1e-15)
    if (norm(mesh1.vertices(bindices1(1),1:2) - mesh2.vertices(bindices2(end),1:2)) > 1e-15)
        error('Meshes are not mergeable!');
    else
        bindices2 = bindices2(end:-1:1);
    end
end

v = mesh1.vertices(:,3:4); 	
v(v == bindex1) = 0;
mesh1.vertices(:,3:4) = v;

mesh2.vertices(bindices2,:) = mesh2.vertices(bindices2,:)*0;

n_boundary_vertices = length(bindices1);
n_boundaries_2 = size(mesh2.boundaries,2);
nvertices1 = size(mesh1.vertices,1);
nvertices2 = size(mesh2.vertices,1);
mesh2.elements(:,1:3) = mesh2.elements(:,1:3) + nvertices1;
elements2 = mesh2.elements(:,1:3);
elements2 = elements2(:);

for j = 1:n_boundaries_2
    mesh2.boundaries{j} = mesh2.boundaries{j} + nvertices1;
end

elements_containing_vertex = cell(nvertices1 + nvertices2,1);

for i = 1:nvertices1
    elements_containing_vertex{i} = mesh1.elements_containing_vertex{i};
end

for i = nvertices1+1:nvertices1+nvertices2
    elements_containing_vertex{i} = mesh2.elements_containing_vertex{i-nvertices1};
end

for i = 1:n_boundary_vertices
    elements_containing_vertex{i} = [mesh1.elements_containing_vertex{bindices1(i)};
                                     mesh2.elements_containing_vertex{bindices2(i)}];
    elements2(elements2 == bindices2(i) + nvertices1) = bindices1(i);
    for j = 1:n_boundaries_2
        mesh2.boundaries{j}(mesh2.boundaries{j} == bindices2(i) + nvertices1) = bindices1(i);
    end
end

n1 = size(mesh2.elements,1);
mesh2.elements(:,1:3) = reshape(elements2,n1,3);

boundaries = cell(size(boundary_map,1),1);

vertices_boundary_1 = mesh1.vertices(:,3:4);
vertices_boundary_2 = mesh2.vertices(:,3:4);
vertices_boundary_1 = -vertices_boundary_1(:);
vertices_boundary_2 = -vertices_boundary_2(:);
% attention: the order might not be consistent between the two meshes
for i = 1:size(boundary_map,1)
    b = [];
    if (boundary_map(i,1) ~= 0)
        b = [b;mesh1.boundaries{boundary_map(i,1)}];
        el1 = mesh1.elements(:,end) == boundary_map(i,1);
        mesh1.elements(el1,end) = i;
        vertices_boundary_1(vertices_boundary_1 == -boundary_map(i,1)) = i;
    end
    if (boundary_map(i,2) ~= 0)
        b = [b;mesh2.boundaries{boundary_map(i,2)}];
        el2 = mesh2.elements(:,end) == boundary_map(i,2);
        mesh2.elements(el2,end) = i;
        vertices_boundary_2(vertices_boundary_2 == -boundary_map(i,2)) = i;
    end
    boundaries{i} = b;
end

mesh1.vertices(:,3:4) = abs(reshape(vertices_boundary_1,nvertices1,2));
mesh2.vertices(:,3:4) = abs(reshape(vertices_boundary_2,nvertices2,2));

mesh.vertices = [mesh1.vertices;mesh2.vertices];
mesh.elements = [mesh1.elements;mesh2.elements];
mesh.elements_containing_vertex = elements_containing_vertex;
mesh.boundaries = boundaries;

% remove zero vertices and scale other indices
bindices2 = unique(bindices2);
for i = 1:length(bindices2)
   newindex = bindices2(i) + nvertices1 - i + 1;
   mesh.vertices(newindex,:) = [];
   aux = mesh.elements(:,1:3);
   el = aux > newindex;
   aux(el) = aux(el) - 1;
   mesh.elements(:,1:3) = aux;
   mesh.elements_containing_vertex{newindex} = [];
   for j = 1:size(boundary_map,1)
       aux = mesh.boundaries{j};
       el = aux > newindex;
       aux(el) = aux(el) - 1;
       mesh.boundaries{j} = aux;
   end
end

mesh.xp = min(mesh1.xp,mesh2.xp);
mesh.yp = min(mesh1.yp,mesh2.yp);
mesh.h = max(mesh1.h,mesh2.h);
mesh.L = max(mesh.vertices(:,1)) - mesh.xp;
mesh.H = max(mesh.vertices(:,2)) - mesh.yp;
mesh.triang = triangulation(mesh.elements(:,1:3),mesh.vertices(:,1),mesh.vertices(:,2));
mesh.type = 'unstructured';
