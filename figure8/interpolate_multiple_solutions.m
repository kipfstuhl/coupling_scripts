function [interpsol] = interpolate_multiple_solutions(solutions,fespace_u,fespace_p)

nsol = size(solutions,2);

mesh = fespace_u.mesh;

X = mesh.X;
Y = mesh.Y;

nx = size(X,1);
ny = size(X,2);

u1 = zeros(nx,ny);
u2 = zeros(nx,ny);
p = zeros(nx,ny);

for i = 1:nx
    for j = 1:ny
        x = X(i,j);
        y = Y(i,j);
        for k = 1:nsol
            curmesh = solutions{k}.fespace_u.mesh;
            if (x >= curmesh.xp && x <= (curmesh.xp + curmesh.L) && ...
                y >= curmesh.yp && y <= (curmesh.yp + curmesh.H))
                u1(i,j) = evaluate_fe_function(solutions{k}.u1,solutions{k}.fespace_u,[x;y]);
                u2(i,j) = evaluate_fe_function(solutions{k}.u2,solutions{k}.fespace_u,[x;y]);
                p(i,j) = evaluate_fe_function(solutions{k}.p,solutions{k}.fespace_p,[x;y]);
            end
        end
    end
end

interpsol.u1 = u1;
interpsol.u2 = u2;
interpsol.p = p;
interpsol.fespace_u = fespace_u;
interpsol.fespace_p = fespace_p;
