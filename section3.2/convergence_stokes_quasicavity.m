clear all
clc

u1ex = @(x) sin(x(1,:)*pi).*x(2,:);
u2ex = @(x) -pi*0.5*cos(x(1,:)*pi).*x(2,:).^2;
pex = @(x) -0.5*x(1,:)^2;

u1exdx = @(x) pi*cos(x(1,:)*pi).*x(2,:);
u1exdy = @(x) sin(x(1,:)*pi);
u1exdxdx = @(x) -pi^2*sin(x(1,:)*pi).*x(2,:);
u1exdydy = @(x) 0*x(1,:);

u2exdx = @(x) 0.5*pi^2*sin(x(1,:)*pi).*x(2,:).^2;
u2exdy = @(x) -pi*cos(x(1,:)*pi).*x(2,:);
u2exdxdx = @(x) 0.5*pi^3*cos(x(1,:)*pi).*x(2,:).^2;
u2exdydy = @(x) -pi*cos(x(1,:)*pi);

graduex = @(x) [u1exdx(x) u1exdy(x);u2exdx(x) u2exdy(x)];

pexdx = @(x) -x(1,:);
pexdy = @(x) 0;

mu = 1;
f = @(x) [-mu*(u1exdxdx(x)+u1exdydy(x)) + pexdx(x);
    -mu*(u2exdxdx(x)+u2exdydy(x)) + pexdy(x)];

dirichlet_functions = @(x) [u1ex(x) u2ex(x);
    u1ex(x)*0 u2ex(x)*0;
    u1ex(x) u2ex(x);
    u1ex(x) u2ex(x)]';
neumann_functions = @(x) [mu*graduex(x)*[0;-1]-pex(x)*[0;-1], ...
    mu*graduex(x)*[1;0]-pex(x)*[1;0], ...
    mu*graduex(x)*[0;1]-pex(x)*[0;1], ...
    mu*graduex(x)*[-1;0]-pex(x)*[-1;0]];

L = 1;
H = 1;

n_elements_x = 200;
n_elements_y = 200;

mesh = create_mesh(0,0,1,1,n_elements_x,n_elements_y);
% 
bc_flags = [1 1 1 1];

fespace_u = create_fespace(mesh,'P2',bc_flags);
fespace_p = create_fespace(mesh,'P1',bc_flags);

dirichlet_functions = @(x) [0 0;
    0 0;
    0.5*(1-cos(2*pi*x(1))) 0;
    0 0]';

neumann_functions = @(x) [0 0;
    0 0;
    0 0;
    0 0]';

f = @(x) [x(1,:)*0;x(1,:)*0];

[A,b] = assembler_steady_stokes(fespace_u,fespace_p,f,mu,dirichlet_functions,...
    neumann_functions);

exactsol = solve_fluid_system(A,b,fespace_u,fespace_p);

u1ex = @(x) evaluate_fe_function(exactsol.u1,fespace_u,x);
u2ex = @(x) evaluate_fe_function(exactsol.u2,fespace_u,x);
pex = @(x) evaluate_fe_function(exactsol.p,fespace_p,x);

graduex = @(x) [evaluate_fe_function_gradient(exactsol.u1,fespace_u,x)';
                evaluate_fe_function_gradient(exactsol.u2,fespace_u,x)'];

nelements = [16 32 64 128];

err = [];
h = [];
for i = nelements
    mesh = create_mesh(0,0,1,1,i,i);
    h = [h mesh.h];
    fespace_u = create_fespace(mesh,'P2',bc_flags);
    fespace_p = create_fespace(mesh,'P1',bc_flags);

    [A,b] = assembler_steady_stokes(fespace_u,fespace_p,f,mu,dirichlet_functions,...
        neumann_functions);


    sol = solve_fluid_system(A,b,fespace_u,fespace_p);
    
    err = [err; compute_H1_error_velocity(fespace_u,sol,@(x) [u1ex(x);u2ex(x)], @(x) graduex(x))]
    % err = [err; compute_L2_error_velocity(fespace_u,sol,@(x) [u1ex(x);u2ex(x)])];
end

loglog(h,err*0.1)
hold on
loglog(h,h.^2)