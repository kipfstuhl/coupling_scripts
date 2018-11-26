function [c,expc] = plot_error_bifurcation(ref,freqq,varargin)
meshes = cell(3,1);

mesh_labels = {'inflow_distorted.msh', ...
    'outflow1_distorted.msh', ...
    'outflow2_distorted.msh'};

% load meshes
for i = 1:3
    meshes{i} = read_mesh(['../meshes/refinement',num2str(ref),'/',mesh_labels{i}]);
end
fespace_us = cell(3,1);
fespace_ps = cell(3,1);

bc_matrix = [1 0 0 1 1; 1 0 1 0 0; 0 1 0 1 0];

% build fespaces
for i = 1:3
    fespace_us{i} = create_fespace(meshes{i},'P2',bc_matrix(i,:));
    fespace_ps{i} = create_fespace(meshes{i},'P1',[]);
end

r_in = 0.5;
U = 1;
dir0 = @(x) zeros(2,6);
dir_in = @(x) [0 0; 0 0; 0 0;0 0; (r_in^2 - x(2)^2)/r_in^2 * U 0]';
neu0 = dir0;


As = cell(3,1);
bs = cell(3,1);

[As{1},bs{1}] = assembler_steady_navier_stokes(fespace_us{1},fespace_ps{1},[0;0],1,dir_in,neu0);
[As{2},bs{2}] = assembler_steady_navier_stokes(fespace_us{2},fespace_ps{2},[0;0],1,dir0,neu0);
[As{3},bs{3}] = assembler_steady_navier_stokes(fespace_us{3},fespace_ps{3},[0;0],1,dir0,neu0);

% find interface 1
b1 = meshes{1}.boundaries{2};
x1 = meshes{1}.vertices(b1(end,end),1:2);
x2 = meshes{1}.vertices(b1(1,1),1:2);
dif = x2 - x1;
l1 = norm(dif);
t1 = dif/l1;
n1 = [t1(2) -t1(1)];

% find interface 2
b2 = meshes{2}.boundaries{4};
x1 = meshes{2}.vertices(b2(1,1),1:2);
x2 = meshes{2}.vertices(b2(end,end),1:2);
dif = x2 - x1;
l2 = norm(dif);
t2= dif/l2;
n2 = [t2(2) -t2(1)];

% find interface 3
b3 = meshes{3}.boundaries{5};
x1 = meshes{3}.vertices(b3(1,1),1:2);
x2 = meshes{3}.vertices(b3(end,end),1:2);
dif = x2 - x1;
l3 = norm(dif);
t3= dif/l3;
n3 = [t3(2) -t3(1)];

ts = [t1;t2;t3];
ns = [n1;n2;n3];
ls = [l1 l2 l3];

nfreqs = [1 1 1] * freqq;

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
    -3 0 5];

n_boundaries = [5 5 5];

xc = x1;

B = [];

prods = zeros(3);

for j = 1:3 % interface
    for i = 0:nfreqs(j)
        if (i == 0)
            bx = zeros(totalnodes,1);
            by = zeros(totalnodes,1);
            f = @(x,xp,l) 1;
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
        else
            bx = zeros(totalnodes,2);
            by = zeros(totalnodes,2);
            f = @(x,xp,l) sin(i*sqrt((x(1,:)-xp(1)).^2+(x(2,:)-xp(2)).^2)/l*pi);
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
            
            f = @(x,xp,l) cos(i*sqrt((x(1,:)-xp(1)).^2+(x(2,:)-xp(2)).^2)/l*pi);
            for k = 1:3 % domains
                flags = zeros(1,n_boundaries(k));
                si = sign(connectivity(j,k));
                if (si ~= 0)
                    flags(abs(connectivity(j,k))) = si;
                    b1 = zeros(n_nodes_us{k},1);
                    b2 = zeros(n_nodes_us{k},1);
                    b1 = apply_neumann_bc(b1,fespace_us{k},@(x) flags * f(x,xc,ls(j)),8);
                    bx(indices{k},2) = bx(indices{k},2) + [b1;b2];
                    by(indices{k},2) = by(indices{k},2) + [b2;b1];
                end
            end
        end
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
    cur = cur + n_nodes_tot{i};
end
B(dir_indices,:) = 0;

nlamb = size(B,2);

jac_block11 = @(u) [];
A = @(u) [];
for i = 1:3
    A = @(u) blkdiag(A(u),As{i}(u(indices{i})));
    jac_block11 = @(u) blkdiag(jac_block11(u),build_jac_navier_stokes(As{i},u(indices{i}),fespace_us{i}));
end

x0 = zeros(length(bs{1}) + length(bs{2}) + length(bs{3}),1);
x0 = [x0;zeros(nlamb,1)];

jac = @(u) [jac_block11(u) B; Bt sparse(nlamb,nlamb)];

H =@(u) [A(u) B;Bt sparse(nlamb,nlamb)];

rhs = [bs{1};bs{2};bs{3};zeros(nlamb,1)];

% solve system with newton's method
ff = @(u) H(u)*u-rhs;
tol = 1e-8;
maxit = 20;

[sol,er,it] = solve_with_newtons_method(ff,rhs,jac,tol,maxit);
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

load(['solutions/u1s_ref',num2str(ref),'.mat']);
load(['solutions/u2s_ref',num2str(ref),'.mat']);
load(['solutions/ps_ref',num2str(ref),'.mat']);

meshes_fine = cell(3,1);

ref = 6;
% now we interpolate on finer meshes for better plot resolution
% load fine meshes
for i = 1:3
    meshes_fine{i} = read_mesh(['../meshes/refinement',num2str(ref),'/',mesh_labels{i}]);
end

fespace_fine_us = cell(3,1);

bc_matrix = [1 0 0 1 1; 1 0 1 0 0; 0 1 0 1 0];

% % create fine fespaces
for i = 1:3
    fespace_fine_us{i} = create_fespace(meshes_fine{i},'P2',bc_matrix(i,:));
end

% norm of exact solution
% normu = 8.443059661260090;

% bcmatrix = [1 0 0 1 1; 1 1 1 0 0; 0 1 1 1 0];
minlog = Inf;
maxlog = -Inf;
Ns = cell(3,1);
for i = 1:3
    indu = find_dirichlet_indices(fespace_us{i});
        
    dif1 = abs(sols{i}.u1 - u1s{i});    
    dif2 = abs(sols{i}.u2 - u2s{i});
    
    N = sqrt(dif1.^2 + dif2.^2);
    N(indu) = 1e-16;

    Ns{i} = log(N);

    minlog = min(min(Ns{i}),minlog);
    maxlog = max(max(Ns{i}),maxlog);
end

if (nargin > 3)
    minlog = varargin{1};
    maxlog = varargin{2};
end

interpsol = cell(3,1);
c = linspace(minlog,maxlog,30);

% interpolate on fine fespaces
for i = 1:3
    interpsol{i} = interp_on_fespace(fespace_us{i},Ns{i},fespace_fine_us{i});
%     [~,h] = tricontf(fespace_fine_us{i}.mesh.vertices(:,1),fespace_fine_us{i}.mesh.vertices(:,2), ...
%        fespace_fine_us{i}.mesh.elements(:,1:3),interpsol{i}(1:size(fespace_fine_us{i}.mesh.vertices(:,1),1)),c);
    h = trisurf(fespace_fine_us{i}.mesh.elements(:,1:3),fespace_fine_us{i}.mesh.vertices(:,1), ...
        fespace_fine_us{i}.mesh.vertices(:,2),interpsol{i}(1:size(fespace_fine_us{i}.mesh.vertices(:,1),1)));
    set(h,'edgecolor','none');
    hold on
end

% interpolate on fine fespaces
% for i = 1:3
% %     [~,h] = tricontf(fespace_fine_us{i}.mesh.vertices(:,1),fespace_fine_us{i}.mesh.vertices(:,2), ...
% %        fespace_fine_us{i}.mesh.elements(:,1:3),interpsol{i}(1:size(fespace_fine_us{i}.mesh.vertices(:,1),1)),c);
%     h = trisurf(fespace_us{i}.mesh.elements(:,1:3),fespace_us{i}.mesh.vertices(:,1), ...
%         fespace_us{i}.mesh.vertices(:,2),Ns{i}(1:size(fespace_us{i}.mesh.vertices(:,1),1)));
%     set(h,'edgecolor','none');
%     hold on
% end

c = linspace(minlog,maxlog,6);
expc = exp(c);
for i = 1:length(c)
    check = 1;
    count = 0;
    while check
        count = count + 1;
        if (floor(expc(i) * 10^count) > 0)
            check = false;
        end
    end
    if (floor(expc(i) * 10^(count+1) - floor(expc(i) * 10^(count))*10) < 5)
        expc(i) = 10^(-count) * floor(expc(i) * 10^count);
    else
        expc(i) = 10^(-count) * (floor(expc(i) * 10^count) + 1);
    end
    
end

colorbar('YTick',c,'YTickLabel',expc,'Location','SouthOutside');

vertices = [meshes{1}.vertices;meshes{2}.vertices;meshes{3}.vertices];

xm = min(vertices(:,1));
xM = max(vertices(:,1));
ym = min(vertices(:,2));
yM = max(vertices(:,2));

axis equal
axis([xm xM ym yM])
set(gca,'color','none')
set(gca,'Visible','off')
hold off

