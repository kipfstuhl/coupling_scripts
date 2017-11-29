% plot first solution

n_nodes_u = size(fespace1_u.nodes,1);

u11 = sol(1:n_nodes_u);
u22 = sol(n_nodes_u+1:2*n_nodes_u);

n_vertices = size(fespace1_u.mesh.vertices,1);
n1 = size(fespace1_u.mesh.X,1);
n2 = size(fespace1_u.mesh.X,2);
U1 = reshape(u11(1:n_vertices),n1,n2);
U2 = reshape(u22(1:n_vertices),n1,n2);

N = sqrt(U1.^2+U2.^2);

[~,c1] = contourf(fespace1_u.mesh.X,fespace1_p.mesh.Y,N);
% c.LineStyle = 'none';
% shading interp
% h = colorbar;
% if (nargin >= 5)
%     if (mlim == 0 && Mlim == 0)
%         mlim = 0; 
%         Mlim = 1;
%     end
%     set(h,'YLim',[mlim*0.9 Mlim*1.1])
%     caxis([ mlim*0.9 Mlim*1.1 ])
%     c.LevelList = linspace(mlim,Mlim,20);
% else
%     m = min(min(N));
%     M = max(max(N));
%     c.LevelList = linspace(m,M,20);
% end
% hold on