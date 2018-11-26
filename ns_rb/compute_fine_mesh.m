clear all
close all
clc

mesh1 = read_mesh('../meshes/inflow_distorted_conforming.msh');
mesh2 = read_mesh('../meshes/outflow1_distorted_conforming.msh');
mesh3 = read_mesh('../meshes/outflow2_distorted_conforming.msh');
%%
boundary_map = [1 0; 2 0; 3 2; 0 3; 0 4; 5 5];
m1 = merge_meshes(mesh2,mesh3,4,1,boundary_map);

disp('Done 1')
boundary_map = [1 1; 0 2; 0 3; 0 4; 3 5; 4 0];
m2 = merge_meshes(mesh1,m1,2,6,boundary_map);

disp('Done 2')
% boundary_map = [1 1; 2 0; 3 0; 4 0; 5 3; 0 4];
% m2 = merge_meshes(m1,mesh1,6,6,boundary_map);


% draw_mesh(m2)
% hold on
% n_vertices = size(m2.vertices,1);
% index = 6;
% for i = 1:n_vertices
%     vv = m2.vertices(i,:);
%     if (vv(3) == index || vv(4) == index)
%        plot(vv(1),vv(2),'.r','Markersize',10) 
%     end
% end

fine_mesh = m2;
save('solutions/fine_mesh.mat','fine_mesh');