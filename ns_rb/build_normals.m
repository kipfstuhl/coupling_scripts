function mesh = build_normals(mesh)

nboundaries = size(mesh.boundaries,2);
normals = cell(nboundaries,1);

for i = 1:nboundaries
    curboundary = mesh.boundaries{i};
    nelements = size(curboundary,1);
    norml = zeros(nelements,2);
    
    for j = 1:nelements
        x1 = mesh.vertices(curboundary(j,1),1:2);
        x2 = mesh.vertices(curboundary(j,2),1:2);
        tangent = (x1-x2)/norm(x1-x2);
        norml(j,:) = [tangent(2) -tangent(1)];
    end
    normals{i} = norml;
end

mesh.normals = normals;