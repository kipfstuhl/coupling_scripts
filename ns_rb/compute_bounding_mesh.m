function [bounding_mesh] = compute_bounding_mesh(mesh,nx,ny)

xp = mesh.xp;
yp = mesh.yp;
L = mesh.L;
H = mesh.H;

bounding_mesh = create_mesh(xp,yp,L,H,nx,ny);
n_elements_coarse = size(bounding_mesh.elements,1);

map_elements = cell(n_elements_coarse,1);

vertices = mesh.vertices;

n_vertices = size(vertices,1);

for i = 1:n_vertices
    index = find_element_containing_point(bounding_mesh,mesh.vertices(i,1:2)');
    map_elements{index} = [map_elements{index};i];
end

for i = 1:n_elements_coarse
    els = bounding_mesh.elements(i,1:3);
    for j = 1:3
        index = find_element_containing_point(mesh,bounding_mesh.vertices(els(j),1:2)');
        if (index ~= -1)
            map_elements{i} = [map_elements{i}; ...
                mesh.elements(index,1); ...
                mesh.elements(index,2); ...
                mesh.elements(index,3)];
        end
    end
end

for i = 1:n_elements_coarse
    map_elements{i} = unique(map_elements{i});
end

bounding_mesh.map_elements = map_elements;

