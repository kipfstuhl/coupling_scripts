clear all
close all
clc

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex')

hold on
mesh_in = read_mesh('../meshes/inflow_symmetric.msh');
% mesh_in = build_normals(mesh_in);
draw_mesh(mesh_in);
mesh_out1 = read_mesh('../meshes/outflow1_symmetric.msh');
% mesh_out1 = build_normals(mesh_out1);
draw_mesh(mesh_out1);
mesh_out2 = read_mesh('../meshes/outflow2_symmetric.msh');
% mesh_out2 = build_normals(mesh_out2);
meshes = {mesh_in,mesh_out1,mesh_out2};
draw_mesh(mesh_out2);

vertices = [mesh_in.vertices;mesh_out1.vertices;mesh_out2.vertices];

xm = min(vertices(:,1));
xM = max(vertices(:,1));
ym = min(vertices(:,2));
yM = max(vertices(:,2));

axis([xm xM ym yM])

fespace_us = cell(3,1);
fespace_ps = cell(3,1);

% fespaces inlet
fespace_us{1} = create_fespace(mesh_in,'P2',[1 0 0 1 1]);

% fespaces outlet1
fespace_us{2} = create_fespace(mesh_out1,'P2',[1 1 1 0 0]);

% fespaces outlet2
fespace_us{3} = create_fespace(mesh_out2,'P2',[0 1 1 1 0]);

ff = @(x) x(2).*sin(x(1));
ff = @(x) (x(2).^2 -2*x(1).^2)*0.1;

U = 1;
dir0 = @(x) zeros(1,6);
dir_in = @(x) ones(1,6)*ff(x);
neu0 = dir0;

As = cell(3,1);
bs = cell(3,1);

[As{1},bs{1}] = assembler_poisson(fespace_us{1},[0;0],1,dir_in,neu0);
[As{2},bs{2}] = assembler_poisson(fespace_us{2},[0;0],1,dir_in,neu0);
[As{3},bs{3}] = assembler_poisson(fespace_us{3},[0;0],1,dir_in,neu0);

nfreq = 5;
nfourier = nfreq*2+1;

% find interface 1
b1 = mesh_in.boundaries{2};
x1 = mesh_in.vertices(b1(1,1),1:2);
x2 = mesh_in.vertices(b1(end,end),1:2);
dif = x2 - x1;
l1 = norm(dif);
t1 = dif/l1;
n1 = [t1(2) -t1(1)];

x1
x2

% find interface 2
b2 = mesh_out1.boundaries{4};
x2 = mesh_out1.vertices(b2(1,1),1:2);
x1 = mesh_out1.vertices(b2(end,end),1:2);
dif = x2 - x1;
l2 = norm(dif);
t2= dif/l2;
n2 = [t2(2) -t2(1)];
x1
x2

plot(x1(1),x1(2),'.','Markersize',20)
plot(x2(1),x2(2),'.','Markersize',20)


% find interface 3
b3 = mesh_out2.boundaries{5};
x1 = mesh_out2.vertices(b3(1,1),1:2);
x2 = mesh_out2.vertices(b3(end,end),1:2);
dif = x2 - x1;
l3 = norm(dif);
t3= dif/l3;
n3 = [t3(2) -t3(1)];
x1
x2

ts = [t1;t2;t3];
ns = [n1;n2;n3];

ls = [l1 l2 l3];
n_nodes_us = cell(3,1);
n_nodes_ps = cell(3,1);
n_nodes_tot = cell(3,1);
indices = cell(3,1);
totalnodes = 0;
cur = 0;
for i = 1:3
    n_nodes_us{i} = size(fespace_us{i}.nodes,1);
    indices{i} = cur+1:cur+n_nodes_us{i};
    cur = cur + n_nodes_us{i};
    totalnodes = totalnodes + n_nodes_us{i};
end

connectivity = [-2 5 0;
    0 4 -1;
    3 0 -5];

n_boundaries = [5 5 5];

xc = x1;

B = [];

% add "constants"
for i = 1:3 % interface 1
    b = zeros(totalnodes,1);
    for j = 1:3 % interface 2
        prod = abs(ts(i,:)*ns(j,:)');
        for k = 1:3 % domains
            flags = zeros(1,n_boundaries(k));
            si = sign(connectivity(j,k));
            if (si ~= 0)
                flags(abs(connectivity(j,k))) = si;
                b1 = zeros(n_nodes_us{k},1);
                b1 = apply_neumann_bc(b1,fespace_us{k},@(x) prod*flags,8);
                b(indices{k}) = b(indices{k}) + b1;
                %                 if (j == 2 && k == 3)
                %                     keyboard
                %                 end
  
                bbb = meshes{k}.boundaries{abs(connectivity(j,k))};
                bbb = [bbb(:,1);bbb(end,end)];
                xxx = meshes{k}.vertices(bbb,1:2);
                
                plot3(xxx(:,1),xxx(:,2),prod*si*xxx(:,1).^0,'Linewidth',1);
                hold on
                pause()
            end
        end
    end
    B = [B b];
end

for i = 1:nfreq
    for j = 1:3 % interface
        b = zeros(totalnodes,1);
        f = @(x,xp,l) sin(i*sqrt((x(1,:)-xp(1)).^2+(x(2,:)-xp(2)).^2)/l*pi/2);
        for k = 1:3 % domains
            flags = zeros(1,n_boundaries(k));
            si = sign(connectivity(j,k));
            if (si ~= 0)
                flags(abs(connectivity(j,k))) = si;
                b1 = zeros(n_nodes_us{k},1);
                b1 = apply_neumann_bc(b1,fespace_us{k},@(x) flags * f(x,xc',ls(j)),8);
                b(indices{k},1) = b(indices{k},1) + b1;
                
                disp(['i = ', num2str(i), ...
                    ' j = ',num2str(j), ...
                    ' k = ',num2str(k)])
                bbb = meshes{k}.boundaries{abs(connectivity(j,k))};
                bbb = [bbb(:,1);bbb(end,end)];
                xxx = meshes{k}.vertices(bbb,1:2);
                
                plot3(xxx(:,1),xxx(:,2),f(xxx(:,1:2)',xc,ls(j)),'Linewidth',1);
                hold on
                pause()
            end  
        end
        %             f = @(x,xp,l) cos(i*norm(x-xp)/l*pi);
        %             for k = 1:3 % domains
        %                 flags = zeros(1,n_boundaries(k));
        %                 si = sign(connectivity(j,k));
        %                 if (si ~= 0)
        %                     flags(abs(connectivity(j,k))) = si;
        %                     b1 = zeros(n_nodes_us{k},1);
        %                     b2 = zeros(n_nodes_us{k},1);
        %                     b1 = apply_neumann_bc(b1,fespace_us{k},@(x) flags * f(x,xc,ls(j)),8);
        %                     bx(indices{k},2) = bx(indices{k},2) + [b1;b2];
        %                     by(indices{k},2) = by(indices{k},2) + [b2;b1];
        %                 end
        %             end
        B = [B b];
    end    
end
B = sparse(B);

Bt = B';
dir_indices = [];
cur = 0;
for i = 1:3
    loc_dir_indices = find_dirichlet_indices(fespace_us{i});
    dir_indices = [dir_indices;
        loc_dir_indices + cur];
    for j = loc_dir_indices(:)'
        plot(fespace_us{i}.nodes(j,1),fespace_us{i}.nodes(j,2),'.r','Markersize',10);
    end
    cur = cur + n_nodes_us{i};
end
B(dir_indices,:) = 0;

nlamb = size(B,2);

A = [As{1} sparse(n_nodes_us{1},n_nodes_us{2}) sparse(n_nodes_us{1},n_nodes_us{3}); ...
    sparse(n_nodes_us{2},n_nodes_us{1}) As{2} sparse(n_nodes_us{2},n_nodes_us{3}); ...
    sparse(n_nodes_us{3},n_nodes_us{1}) sparse(n_nodes_us{3},n_nodes_us{2}) As{3}];

H = [A B;Bt sparse(nlamb,nlamb)];
rhs = [bs{1};bs{2};bs{3};zeros(nlamb,1)];

sol = H\rhs;
sols = cell(3,1);

figure
hold on
axis equal
for i = 1:3
    sols{i} = sol(indices{i});
    plot_fe_function(sols{i},fespace_us{i})
end


export_vtk_scalar(sols{1},fespace_us{1},'sol1')
export_vtk_scalar(sols{2},fespace_us{2},'sol2')
export_vtk_scalar(sols{3},fespace_us{3},'sol3')
