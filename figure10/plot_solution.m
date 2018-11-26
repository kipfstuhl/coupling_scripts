clear all
clc

% refinement of the meshes
ref = 3;

% number of frequencies on the interfaces 
freqq = 4;

sols = solve_system_bifurcation_FEM_FEM(ref,freqq);

% now we interpolate on finer meshes for better plot resolution
meshes_fine = cell(3,1);
ref = 6;

mesh_labels = {'inflow_distorted.msh', ...
    'outflow1_distorted.msh', ...
    'outflow2_distorted.msh'};

% load fine meshes
for i = 1:3
     meshes_fine{i} = read_mesh(['../meshes/refinement',num2str(ref),'/',mesh_labels{i}]);
end

fespace_fine_us = cell(3,1);

% create fine fespaces
bc_matrix = [1 0 0 1 1; 1 0 1 0 0; 0 1 0 1 0];
for i = 1:3
    fespace_fine_us{i} = create_fespace(meshes_fine{i},'P2',bc_matrix(i,:));
end

% interpolate on fine fespaces
for i = 1:3
    sols{i}.u1 = interp_on_fespace(sols{i}.fespace_u,sols{i}.u1,fespace_fine_us{i});
    sols{i}.u2 = interp_on_fespace(sols{i}.fespace_u,sols{i}.u2,fespace_fine_us{i});
    sols{i}.fespace_u = fespace_fine_us{i};
end

% plot solutions
hold on
unorm = cell(3,1);
for i = 1:3
   unorm{i} = sqrt(sols{i}.u1.^2 + sols{i}.u2.^2); 
   unorm{i} = unorm{i}(1:size(fespace_fine_us{i}.mesh.vertices(:,1),1));
end

minnorm = 0;
maxnorm = 1.6;

% set the same range of contourlines in all the subdomains
c = linspace(minnorm,maxnorm,20);

for i = 1:3
    [~,h] = tricontf(fespace_fine_us{i}.mesh.vertices(:,1),fespace_fine_us{i}.mesh.vertices(:,2), ...
       fespace_fine_us{i}.mesh.elements(:,1:3),unorm{i}(1:size(fespace_fine_us{i}.mesh.vertices(:,1),1)),c);
    hold on
end

vertices = [meshes_fine{1}.vertices;meshes_fine{2}.vertices;meshes_fine{3}.vertices];

xm = min(vertices(:,1));
xM = max(vertices(:,1));
ym = min(vertices(:,2));
yM = max(vertices(:,2));

axis equal
axis([xm xM ym yM])
set(gca,'color','none')
set(gca,'Visible','off')
hold off
