clear all
close all
clc

mesh_in = read_mesh('../meshes/refinement1/inflow_distorted.msh');
mesh_out1 = read_mesh('../meshes/refinement1/outflow1_distorted.msh');
mesh_out2 = read_mesh('../meshes/refinement1/outflow2_distorted.msh');

hold on

trisurf(mesh_in.elements(:,1:3),mesh_in.vertices(:,1), ...
    mesh_in.vertices(:,2),zeros(size(mesh_in.vertices,1),1),'facecolor',[255 217 230]/255);
trisurf(mesh_out1.elements(:,1:3),mesh_out1.vertices(:,1), ...
    mesh_out1.vertices(:,2),zeros(size(mesh_out1.vertices,1),1),'facecolor',[178 236 174]/255);
trisurf(mesh_out2.elements(:,1:3),mesh_out2.vertices(:,1), ...
    mesh_out2.vertices(:,2),zeros(size(mesh_out2.vertices,1),1),'facecolor',[215 227 255]/255);

vertices = [mesh_in.vertices;mesh_out1.vertices;mesh_out2.vertices];

xm = min(vertices(:,1));
xM = max(vertices(:,1));
ym = min(vertices(:,2));
yM = max(vertices(:,2));

axis equal

axis([xm xM ym yM])

set(gca,'color','none')
set(gca,'Visible','off')

% find interface 1
b1 = mesh_in.boundaries{2};
x1 = mesh_in.vertices(b1(end,end),1:2);
x2 = mesh_in.vertices(b1(1,1),1:2);
plot([x1(1) x2(1)],[x1(2) x2(2)],'Linewidth',3,'Color',[0 167 0]/255);

% find interface 2
b2 = mesh_out1.boundaries{4};
x1 = mesh_out1.vertices(b2(1,1),1:2);
x2 = mesh_out1.vertices(b2(end,end),1:2);
plot([x1(1) x2(1)],[x1(2) x2(2)],'Linewidth',3,'Color',[78 0 255]/255);

% find interface 3
b3 = mesh_out2.boundaries{5};
x1 = mesh_out2.vertices(b3(1,1),1:2);
x2 = mesh_out2.vertices(b3(end,end),1:2);
plot([x1(1) x2(1)],[x1(2) x2(2)],'Linewidth',3,'Color',[255 0 0]/255);

hold off


