function [c,expc] = plot_error_bifurcation(ref,freqq,interp_ref,minlog,maxlog)

sols = solve_system_bifurcation_FEM_FEM(ref,freqq);

load(['data_figure11/u1s_ref',num2str(ref),'.mat'],'u1s');
load(['data_figure11/u2s_ref',num2str(ref),'.mat'],'u2s');

meshes_fine = cell(3,1);

% now we interpolate on anoter mesh to have consistent resolution
mesh_labels = {'inflow_distorted.msh', ...
    'outflow1_distorted.msh', ...
    'outflow2_distorted.msh'};
for i = 1:3
    meshes_fine{i} = read_mesh(['../meshes/refinement',num2str(interp_ref),'/',mesh_labels{i}]);
end

fespace_fine_us = cell(3,1);

bc_matrix = [1 0 0 1 1; 1 0 1 0 0; 0 1 0 1 0];
% create fine fespaces
for i = 1:3
    fespace_fine_us{i} = create_fespace(meshes_fine{i},'P2',bc_matrix(i,:));
end

Ns = cell(3,1);
for i = 1:3
    indu = find_dirichlet_indices(sols{i}.fespace_u);
        
    dif1 = abs(sols{i}.u1 - u1s{i});    
    dif2 = abs(sols{i}.u2 - u2s{i});
    
    N = sqrt(dif1.^2 + dif2.^2);
    N(indu) = 1e-16;

    Ns{i} = log(N);
end

interpsol = cell(3,1);

% interpolate on fine fespaces
for i = 1:3
    interpsol{i} = interp_on_fespace(sols{i}.fespace_u,Ns{i},fespace_fine_us{i});
    h = trisurf(fespace_fine_us{i}.mesh.elements(:,1:3),fespace_fine_us{i}.mesh.vertices(:,1), ...
        fespace_fine_us{i}.mesh.vertices(:,2),interpsol{i}(1:size(fespace_fine_us{i}.mesh.vertices(:,1),1)));
    set(h,'edgecolor','none');
    hold on
end

% create label for the colorbar (we do it manually since the colorbar is
% logarithmic)
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

vertices = [meshes_fine{1}.vertices; ...
            meshes_fine{2}.vertices; ...
            meshes_fine{3}.vertices];

xm = min(vertices(:,1));
xM = max(vertices(:,1));
ym = min(vertices(:,2));
yM = max(vertices(:,2));

axis equal
axis([xm xM ym yM])
set(gca,'color','none')
set(gca,'Visible','off')
hold off

shading interp
colormap jet
view(0,90)
caxis([minlog maxlog])

