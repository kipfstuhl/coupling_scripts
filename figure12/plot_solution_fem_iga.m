% subtract 10 to the norm so that the mesh (plotted at z=0) is visible
h = trisurf(sol_fem.fespace_u.mesh.elements(:,1:3),sol_fem.fespace_u.mesh.vertices(:,1), ...
    sol_fem.fespace_u.mesh.vertices(:,2),sqrt(sol_fem.u1(1:size(sol_fem.fespace_u.mesh.vertices(:,2),1))).^2 + sol_fem.u2(1:size(sol_fem.fespace_u.mesh.vertices(:,2),1)).^2-10);
set(h,'edgecolor','none');
hold on

for i = 1:2
    [eu,ef] = sp_eval(is{i}.vel,is{i}.space_v,is{i}.geometry,[101 101]);
    X = squeeze(ef(1,:,:));
    Y = squeeze(ef(2,:,:));
    Z = squeeze(sqrt(eu(1,:,:).^2 + eu(2,:,:).^2))-10;
    surf(X,Y,Z);
    view(0,90)
end
shading interp

axis equal
set(gca,'color','none')
set(gca,'Visible','off')

% draw meshes
draw_mesh(sol_fem.fespace_u.mesh,[0 0 0])
hold on

vtk_pts = {linspace(0, 1, 51), linspace(0, 1, 51)};
for i = 1:2
    hmsh = hierarchical_mesh(is{i}.msh,[2 2]);
    hmsh_plot_cells(hmsh);   
end

create_inflow_nurb
for i = 1:ni
    for j = 1:nj
        p = coefs(:,i,j);
        plot(p(1),p(2),'.r','Markersize',20);
    end
end

create_outflow2_nurb
for i = 1:ni
    for j = 1:nj
        p = coefs(:,i,j);
        plot(p(1),p(2),'.r','Markersize',20);
    end
end

xm = 0;
xM = 4.16;
ym = -0.6;
yM = 1.15;

axis([xm xM ym yM])