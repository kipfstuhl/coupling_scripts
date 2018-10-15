clear all
close all
clc

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex')

hold on
mesh_in = read_mesh('../meshes/inflow.msh');
% mesh_in = build_normals(mesh_in);
draw_mesh(mesh_in);
mesh_out1 = read_mesh('../meshes/outflow1.msh');
% mesh_out1 = build_normals(mesh_out1);
draw_mesh(mesh_out1);
mesh_out2 = read_mesh('../meshes/outflow2.msh');
% mesh_out2 = build_normals(mesh_out2);
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
fespace_ps{1} = create_fespace(mesh_in,'P1',[]);

% fespaces outlet1
fespace_us{2} = create_fespace(mesh_out1,'P2',[1 0 1 0 0]);
fespace_ps{2} = create_fespace(mesh_out1,'P1',[]);

% fespaces outlet2
fespace_us{3} = create_fespace(mesh_out2,'P2',[0 1 0 1 0]);
fespace_ps{3} = create_fespace(mesh_out2,'P1',[]);

r_in = 0.5;
U = 1;
dir0 = @(x) zeros(2,6);
dir_in = @(x) [0 0; 0 0; 0 0;0 0; (r_in^2 - x(2)^2)/r_in^2 * U 0]';
neu0 = dir0;

As = cell(3,1);
bs = cell(3,1);

[As{1},bs{1}] = assembler_steady_stokes(fespace_us{1},fespace_ps{1},[0;0],1,dir_in,neu0);
[As{2},bs{2}] = assembler_steady_stokes(fespace_us{2},fespace_ps{2},[0;0],1,dir0,neu0);
[As{3},bs{3}] = assembler_steady_stokes(fespace_us{3},fespace_ps{3},[0;0],1,dir0,neu0);

nfreq = 14;
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
    n_nodes_ps{i} = size(fespace_ps{i}.nodes,1);
    n_nodes_tot{i} = 2*n_nodes_us{i} + n_nodes_ps{i};
    indices{i} = cur+1:cur+2*n_nodes_us{i};
    cur = cur + n_nodes_tot{i};
    totalnodes = totalnodes + n_nodes_tot{i};
end

connectivity = [2 -5 0;
    0 4 -1;
    3 0 -5];

n_boundaries = [5 5 5];

xc = x1;

B = [];

% add "constants"
% for i = 1:3 % interface 1
%     bx = zeros(totalnodes,1);
%     by = zeros(totalnodes,1);
%     for j = 1:3 % interface 2
%         prod = abs(ts(i,:)*ns(j,:)');
%         disp(['i = ', num2str(i),' j = ',num2str(j),' prod = ',num2str(prod)])
%         for k = 1:3 % domains
%             flags = zeros(1,n_boundaries(k));
%             si = sign(connectivity(j,k));
%             if (si ~= 0)
%                 flags(abs(connectivity(j,k))) = si;
%                 b1 = zeros(n_nodes_us{k},1);
%                 b2 = zeros(n_nodes_us{k},1);
%                 b1 = apply_neumann_bc(b1,fespace_us{k},@(x) prod*flags,8);
%                 bx(indices{k}) = bx(indices{k}) + [b1;b2];
%                 by(indices{k}) = by(indices{k}) + [b2;b1];
%             end
%         end
%     end
%     B = [B bx by];
% end

for i = 0:nfreq
    for j = 1:3 % interface
        bx = zeros(totalnodes,1);
        by = zeros(totalnodes,1);
        if (i == 0)
            f = @(x,xp,l) 1;
        else
            f = @(x,xp,l) sin(i*sqrt((x(1,:)-xp(1)).^2+(x(2,:)-xp(2)).^2)/l*pi/2);
        end
        for k = 1:3 % domains
            flags = zeros(1,n_boundaries(k));
            si = sign(connectivity(j,k));
            if (si ~= 0)
                flags(abs(connectivity(j,k))) = si;
                b1 = zeros(n_nodes_us{k},1);
                b2 = zeros(n_nodes_us{k},1);
                b1 = apply_neumann_bc(b1,fespace_us{k},@(x) flags * f(x,xc,ls(j)),8);
                bx(indices{k},1) = bx(indices{k},1) + [b1;b2];
                by(indices{k},1) = by(indices{k},1) + [b2;b1];
            end
        end
        
%         f = @(x,xp,l) cos(i*norm(x-xp)/l*pi);
%         for k = 1:3 % domains
%             flags = zeros(1,n_boundaries(k));
%             si = sign(connectivity(j,k));
%             if (si ~= 0)
%                 flags(abs(connectivity(j,k))) = si;
%                 b1 = zeros(n_nodes_us{k},1);
%                 b2 = zeros(n_nodes_us{k},1);
%                 b1 = apply_neumann_bc(b1,fespace_us{k},@(x) flags * f(x,xc,ls(j)),8);
%                 bx(indices{k},2) = bx(indices{k},2) + [b1;b2];
%                 by(indices{k},2) = by(indices{k},2) + [b2;b1];
%             end
%         end  
        B = [B bx by];
    end
end
B = sparse(B);

Bt = B';
dir_indices = [];
cur = 0;
for i = 1:3
    loc_dir_indices = find_dirichlet_indices(fespace_us{i});
    dir_indices = [dir_indices;
        loc_dir_indices + cur;
        loc_dir_indices + cur + n_nodes_us{i}];
    for j = loc_dir_indices(:)'
        plot(fespace_us{i}.nodes(j,1),fespace_us{i}.nodes(j,2),'.r','Markersize',10);
    end
    cur = cur + n_nodes_tot{i};
end
B(dir_indices,:) = 0;

nlamb = size(B,2);

A = [As{1} sparse(n_nodes_tot{1},n_nodes_tot{2}) sparse(n_nodes_tot{1},n_nodes_tot{3}); ...
    sparse(n_nodes_tot{2},n_nodes_tot{1}) As{2} sparse(n_nodes_tot{2},n_nodes_tot{3}); ...
    sparse(n_nodes_tot{3},n_nodes_tot{1}) sparse(n_nodes_tot{3},n_nodes_tot{2}) As{3}];

H = [A B;Bt sparse(nlamb,nlamb)];
rhs = [bs{1};bs{2};bs{3};zeros(nlamb,1)];

sol = H\rhs;

sols = cell(3,1);
cur = 0;
for i = 1:3
    sols{i}.u1 = sol(cur+1:cur+n_nodes_us{i});
    sols{i}.u2 = sol(cur+n_nodes_us{i}+1:cur+n_nodes_us{i}*2);
    sols{i}.p = sol(cur+2*n_nodes_us{i}+1:cur+n_nodes_tot{i});
    sols{i}.fespace_u = fespace_us{i};
    sols{i}.fespace_p = fespace_ps{i};
    cur = cur + n_nodes_tot{i};
end

% plot_fe_fluid_function(sols{1},'U')
% hold on
% plot_fe_fluid_function(sols{2},'U')
% plot_fe_fluid_function(sols{3},'U')

plot_fe_function(sols{1}.u1,fespace_us{1})
plot_fe_function(sols{2}.u1,fespace_us{2})
plot_fe_function(sols{3}.u1,fespace_us{3})

export_vtk_fluid(sols{1},'sol1')
export_vtk_fluid(sols{2},'sol2')
export_vtk_fluid(sols{3},'sol3')
