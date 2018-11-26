% TEST_CONVERGENCE: test convergence of NavierStokes solver
clc
clear all
close all

clear problem_data
nurb = nrbsquare([0 0],1,1);
problem_data.geo_name = nurb;

% boundary conditions
problem_data.drchlt_sides = [1 2];
problem_data.nmnn_sides = [3 4];

% Physical parameters
mu = 1;
problem_data.viscosity = @(x, y) mu * ones (size (x));

% Force term
u1ex = @(x,y) cos(y*pi);
u2ex = @(x,y) x.*(x-1);
pex = @(x,y)  0*x;

u1exdx = @(x,y) 0*x;
u1exdy = @(x,y) -pi*sin(y*pi);
u1exdxdx = @(x,y) 0*x;
u1exdydy = @(x,y) -pi^2*cos(y*pi);

u2exdx = @(x,y) 2*x -x.^0;
u2exdy = @(x,y) 0*x;
u2exdxdx = @(x,y) 0*x;
u2exdydy = @(x,y) 2*x.^0;

pexdx = @(x,y) x*0;
pexdy = @(x,y) 0*x;

fx = @(x, y) -mu*(u1exdxdx(x,y)+u1exdydy(x,y)) + (u1ex(x,y).*u1exdx(x,y) + u2ex(x,y).*u1exdy(x,y)) + pexdx(x,y);
fy = @(x, y) -mu*(u2exdxdx(x,y)+u2exdydy(x,y)) + (u1ex(x,y).*u2exdx(x,y) + u2ex(x,y).*u2exdy(x,y)) + pexdy(x,y);


problem_data.f  = @(x, y) cat(1, ...
    reshape (fx (x,y), [1, size(x)]), ...
    reshape (fy (x,y), [1, size(x)]));

% Boundary terms
problem_data.h  = @(x, y, iside) cat(1,reshape(u1ex(x,y),[1,size(x)]), ...
                                       reshape(u2ex(x,y),[1,size(x)]));
problem_data.g  = @(x, y, iside) cat(1,zeros([1,size(x)]),zeros([1,size(x)]));

% Exact solution, to compute the errors
problem_data.velex = @(x, y) cat(1, ...
    reshape (u1ex (x,y), [1, size(x)]), ...
    reshape (u2ex (x,y), [1, size(x)]));

problem_data.gradvelex = @test_navier_stokes_convergence_graduex;

problem_data.pressex = @(x, y) 0*x;

for i = 1:4
    clear method_data
    method_data.element_name = 'th';     % Element type for discretization
    method_data.degree       = [ 3  3];  % Degree of the splines (pressure space)
    method_data.regularity   = [ 2  2];  % Regularity of the splines (pressure space)
    method_data.nsub         = [ 3*2^i  3*2^i];  % Number of subdivisions
    method_data.nquad        = [ 5  5];  % Points for the Gaussian quadrature rule
    
    % Penalization parameter for Nitsche's method
    factor = 10;
    method_data.Cpen = factor*(min(method_data.degree)+1);
    
    solver_data.name = 'newton';
    solver_data.tol = 1e-8;
    solver_data.nmax = 20;
    
    [geometry, msh, space_v, vel, space_p, press] = ...
        solve_navier_stokes (problem_data, method_data,solver_data);
    
    error_l2_p(i) = sp_l2_error (space_p, msh, press, problem_data.pressex)
    [error_h1_v(i), error_l2_v(i)] = ...
        sp_h1_error (space_v, msh, vel, problem_data.velex, problem_data.gradvelex)
    
    msh_prc = msh_precompute(msh);
    h(i) = max(msh_prc.element_size);
end

%%

figure
loglog(h,error_h1_v,'.-','Markersize',10,'Linewidth',1);
title('H1-error on velocity');
hold on
loglog(h,h.^3,'--k');
grid on
xlabel('h');
ylabel('error');
legend('error','h^3');

log(error_h1_v(2:end)./error_h1_v(1:end-1))./log(h(2:end)./h(1:end-1))

figure
loglog(h,error_l2_p,'.-','Markersize',10,'Linewidth',1);
title('L2-error on pressure');
hold on
loglog(h,h.^4*1e-2,'--k');
grid on
xlabel('h');
ylabel('error');
legend('error','h^4');

log(error_l2_p(2:end)./error_l2_p(1:end-1))./log(h(2:end)./h(1:end-1))
