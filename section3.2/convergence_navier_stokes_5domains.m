clear all
close all
clc

% author: Luca Pegolotti on 11/12/2017

set = 1;

if (set == 1)
    U = 1;
    u1ex = @(x) U*sin(x(2,:)*pi);
    u2ex = @(x) 0*x(2,:).^0;
    pex = @(x) -pi^2*U*sin(x(2,:)*pi).*x(1,:);
    
    u1exdx = @(x) 0*x(2,:).^0;
    u1exdy = @(x) pi*U*cos(x(2,:)*pi);
    u1exdxdx = @(x) 0*x(2,:).^0;
    u1exdydy = @(x) -pi^2*U*sin(x(2,:)*pi);
    
    u2exdx = @(x) 0*x(2,:).^0;
    u2exdy = @(x) 0*x(2,:).^0;
    u2exdxdx = @(x) 0*x(2,:).^0;
    u2exdydy = @(x) 0*x(2,:).^0;
    
    graduex = @(x) [u1exdx(x) u1exdy(x);u2exdx(x) u2exdy(x)];
    
    pexdx = @(x)  -pi^2*U*sin(x(2,:)*pi);
    pexdy = @(x)  -pi^3*U*cos(x(2,:)*pi).*x(1,:);
    
    mu = 1;
    nu = mu;
    fun = @(x) [-mu*(u1exdxdx(x)+u1exdydy(x)) + u1ex(x).*u1exdx(x) + u2ex(x).*u1exdy(x) + pexdx(x);
        -mu*(u2exdxdx(x)+u2exdydy(x)) + u1ex(x).*u2exdx(x) + u2ex(x).*u2exdy(x) + pexdy(x)];
    
    dirichlet_functions = @(x) [u1ex(x).*(x(2,:)==0) u2ex(x).*(x(2,:)==0);
        u1ex(x).*(x(1,:)==1) u2ex(x).*(x(1,:)==1);
        u1ex(x).*(x(2,:)==1) u2ex(x).*(x(2,:)==1);
        u1ex(x).*(x(1,:)==0) u2ex(x).*(x(1,:)==0)]';
    neumann_functions = @(x) [(mu*graduex(x)*[0;-1]-pex(x)*[0;-1]).*(x(2,:)==0), ...
        (mu*graduex(x).*[1;0]-pex(x)*[1;0]).*(x(1,:)==1), ...
        (mu*graduex(x).*[0;1]-pex(x)*[0;1]).*(x(2,:)==1), ...
        (mu*graduex(x).*[-1;0]-pex(x)*[-1;0]).*(x(1,:)==0)];
elseif (set == 2)
    U = 1;
    u1ex = @(x) 0*x(2,:).^0;
    u2ex = @(x) U*sin(x(1,:)*pi);
    pex = @(x) -pi^2*U*sin(x(1,:)*pi).*x(2,:);
    
    u1exdx = @(x) 0*x(2,:).^0;
    u1exdy = @(x) 0*x(2,:).^0;
    u1exdxdx = @(x) 0*x(2,:).^0;
    u1exdydy = @(x) 0*x(2,:).^0;
    
    u2exdx = @(x) pi*U*cos(x(1,:)*pi);
    u2exdy = @(x) 0*x(2,:).^0;
    u2exdxdx = @(x) -pi^2*U*sin(x(1,:)*pi);
    u2exdydy = @(x) 0*x(2,:).^0;
    
    graduex = @(x) [u1exdx(x) u1exdy(x);u2exdx(x) u2exdy(x)];
    
    pexdx = @(x)  -pi^3*U*cos(x(1,:)*pi).*x(2,:);
    pexdy = @(x)  -pi^2*U*sin(x(1,:)*pi);
    
    mu = 1;
    nu = mu;
    fun = @(x) [-mu*(u1exdxdx(x)+u1exdydy(x)) + u1ex(x).*u1exdx(x) + u2ex(x).*u1exdy(x) + pexdx(x);
        -mu*(u2exdxdx(x)+u2exdydy(x)) + u1ex(x).*u2exdx(x) + u2ex(x).*u2exdy(x) + pexdy(x)];
    
    dirichlet_functions = @(x) [u1ex(x).*(x(2,:)==0) u2ex(x).*(x(2,:)==0);
        u1ex(x).*(x(1,:)==1) u2ex(x).*(x(1,:)==1);
        u1ex(x).*(x(2,:)==1) u2ex(x).*(x(2,:)==1);
        u1ex(x).*(x(1,:)==0) u2ex(x).*(x(1,:)==0)]';
    neumann_functions = @(x) [(mu*graduex(x)*[0;-1]-pex(x)*[0;-1])*(x(2,:)==0), ...
        (mu*graduex(x)*[1;0]-pex(x)*[1;0])*(x(1,:)==1), ...
        (mu*graduex(x)*[0;1]-pex(x)*[0;1])*(x(2,:)==1), ...
        (mu*graduex(x)*[-1;0]-pex(x)*[-1;0])*(x(1,:)==0)];
elseif (set == 3)
    u1ex = @(x) sin(x(2,:)*pi);
    u2ex = @(x) exp(x(1,:));
    pex = @(x) -0.5*x(1,:).^2;
    
    u1exdx = @(x) 0;
    u1exdy = @(x) pi*cos(x(2,:)*pi);
    u1exdxdx = @(x) 0;
    u1exdydy = @(x) -pi^2*sin(x(2,:)*pi);
    
    u2exdx = @(x) exp(x(1,:));
    u2exdy = @(x) 0;
    u2exdxdx = @(x) exp(x(1,:));
    u2exdydy = @(x) 0;
    
    graduex = @(x) [u1exdx(x) u1exdy(x);u2exdx(x) u2exdy(x)];
    
    pexdx = @(x) -x(1,:);
    pexdy = @(x) 0;
    
    mu = 1;
    nu = mu;
    fun = @(x) [-mu*(u1exdxdx(x)+u1exdydy(x)) + u1ex(x).*u1exdx(x) + u2ex(x).*u1exdy(x) + pexdx(x);
        -mu*(u2exdxdx(x)+u2exdydy(x)) + u1ex(x).*u2exdx(x) + u2ex(x).*u2exdy(x) + pexdy(x)];
    
    dirichlet_functions = @(x) [u1ex(x).*(x(2,:)==0) u2ex(x).*(x(2,:)==0);
        u1ex(x).*(x(1,:)==1) u2ex(x).*(x(1,:)==1);
        u1ex(x).*(x(2,:)==1) u2ex(x).*(x(2,:)==1);
        u1ex(x).*(x(1,:)==0) u2ex(x).*(x(1,:)==0)]';
    neumann_functions = @(x) [(mu*graduex(x)*[0;-1]-pex(x)*[0;-1])*(x(2,:)==0), ...
        (mu*graduex(x)*[1;0]-pex(x)*[1;0])*(x(1,:)==1), ...
        (mu*graduex(x)*[0;1]-pex(x)*[0;1])*(x(2,:)==1), ...
        (mu*graduex(x)*[-1;0]-pex(x)*[-1;0])*(x(1,:)==0)];
    
elseif (set == 4)
    U = 1;
    u1ex = @(x) U*sin(x(2,:)*pi);
    u2ex = @(x) U*sin(x(1,:)*pi);
    pex = @(x) -pi^2*U*sin(x(2,:)*pi).*x(1,:)-pi^2*U*sin(x(1,:)*pi).*x(2,:);
    
    u1exdx = @(x) 0*x(2,:).^0;
    u1exdy = @(x) pi*U*cos(x(2,:)*pi);
    u1exdxdx = @(x) 0*x(2,:).^0;
    u1exdydy = @(x) -pi^2*U*sin(x(2,:)*pi);
    
    u2exdx = @(x) pi*U*cos(x(1,:)*pi);
    u2exdy = @(x) 0*x(2,:).^0;
    u2exdxdx = @(x) -pi^2*U*sin(x(1,:)*pi);
    u2exdydy = @(x) 0*x(2,:).^0;
    
    graduex = @(x) [u1exdx(x) u1exdy(x);u2exdx(x) u2exdy(x)];
    
    pexdx = @(x)  -pi^2*U*sin(x(2,:)*pi)-pi^3*U*cos(x(1,:)*pi).*x(2,:);
    pexdy = @(x)  -pi^3*U*cos(x(2,:)*pi).*x(1,:)-pi^2*U*sin(x(1,:)*pi);
    
    mu = 1;
    nu = mu;
    fun = @(x) [-mu*(u1exdxdx(x)+u1exdydy(x)) + u1ex(x).*u1exdx(x) + u2ex(x).*u1exdy(x) + pexdx(x);
                -mu*(u2exdxdx(x)+u2exdydy(x)) + u1ex(x).*u2exdx(x) + u2ex(x).*u2exdy(x) + pexdy(x)];
    
    dirichlet_functions = @(x) [u1ex(x).*(x(2,:)==0) u2ex(x).*(x(2,:)==0);
        u1ex(x).*(x(1,:)==1) u2ex(x).*(x(1,:)==1);
        u1ex(x).*(x(2,:)==1) u2ex(x).*(x(2,:)==1);
        u1ex(x).*(x(1,:)==0) u2ex(x).*(x(1,:)==0)]';
    neumann_functions = @(x) [(mu*graduex(x)*[0;-1]-pex(x)*[0;-1]).*(x(2,:)==0), ...
        (mu*graduex(x)*[1;0]-pex(x)*[1;0]).*(x(1,:)==1), ...
        (mu*graduex(x)*[0;1]-pex(x)*[0;1]).*(x(2,:)==1), ...
        (mu*graduex(x)*[-1;0]-pex(x)*[-1;0]).*(x(1,:)==0)];
end
n_elementsx = [16 32 64 128];%128];

h = 1./n_elementsx;
errH1u = [];
errL2p = [];
err = [];

xline1 = 0.3;
yline1 = 0.2;

xline2 = 0.6;
yline2 = 0.4;

for nx = n_elementsx
    % create the mesh and fespaces for domain 1
    xp1 = 0;
    yp1 = yline2;
    L1 = 1;
    H1 = 1-yline2;
    
    n1x = nx;
    n1y = floor(n1x*H1);
    mesh1 = create_mesh(xp1,yp1,L1,H1,n1x,n1y);
    
    h_coarse = L1/n1x;
    
    bc_flags = [0 0 1 1];
    fespace1_u = create_fespace(mesh1,'P2',bc_flags);
    fespace1_p = create_fespace(mesh1,'P1',bc_flags);
    
    % create the mesh and fespaces for domain 2
    xp2 = xline2;
    yp2 = 0;
    L2 = 1-xline2;
    H2 = yline2;
    
    n2x = nx/2;
    n2y = nx/2;
    mesh2 = create_mesh(xp2,yp2,L2,H2,n2x,n2y);
    
    bc_flags = [1 0 0 0];
    fespace2_u = create_fespace(mesh2,'P2',bc_flags);
    fespace2_p = create_fespace(mesh2,'P1',bc_flags);
    
    % create the mesh and fespaces for domain 3
    xp3 = xline1;
    yp3 = 0;
    L3 = 1-xline1-L2;
    H3 = yline2;
    
    n3x = nx/2;
    n3y = floor(nx/2*0.9);
    mesh3 = create_mesh(xp3,yp3,L3,H3,n3x,n3y);
    
    bc_flags = [1 0 0 0];
    fespace3_u = create_fespace(mesh3,'P2',bc_flags);
    fespace3_p = create_fespace(mesh3,'P1',bc_flags);
    
    % create the mesh and fespaces for domain 4
    xp4 = 0;
    yp4 = 0;
    L4 = xline1;
    H4 = yline1;
    
    n4x = floor(nx/2*0.9);
    n4y = nx/2;
    mesh4 = create_mesh(xp4,yp4,L4,H4,n4x,n4y);
    
    bc_flags = [1 0 0 1];
    fespace4_u = create_fespace(mesh4,'P2',bc_flags);
    fespace4_p = create_fespace(mesh4,'P1',bc_flags);
    
    % create the mesh and fespaces for domain 5
    xp5 = 0;
    yp5 = yline1;
    L5 = xline1;
    H5 = 1-yline1-H1;
    
    n5x = nx/2;
    n5y = nx/2;
    mesh5 = create_mesh(xp5,yp5,L5,H5,n5x,n5y);
    
    bc_flags = [0 0 0 1];
    fespace5_u = create_fespace(mesh5,'P2',bc_flags);
    fespace5_p = create_fespace(mesh5,'P1',bc_flags);
    
%     draw_multimesh({mesh1, mesh2, mesh3, mesh4, mesh5})
%     pause();

    fespaces_u = {fespace1_u,fespace2_u,fespace3_u,fespace4_u,fespace5_u};
    fespaces_p = {fespace1_p,fespace2_p,fespace3_p,fespace4_p,fespace5_p};

    domain_connectivity = [1 3 3 0 3; 0 4 2 0 0; 0 0 4 2 2; 0 0 0 3 1];
    normals = [-1 0 0 0; 0 0 1 1; 0 -1 1 -1; 0 1 -1 0; 1 1 1 0];
    gausspoints = 4;
    typebasisfunctions = 'polynomial';
    nbasisfunctions = [6 6 6 6];
    [mat,rhs,jac,nsys,nus,nps,indices] = build_coupled_system_navier_stokes(fespaces_u,fespaces_p,fun,nu,dirichlet_functions,neumann_functions,domain_connectivity,normals,nbasisfunctions,gausspoints,typebasisfunctions);
    
    % solve system with newton's method
    f = @(u) mat(u)*u-rhs;
    x0 = zeros(nsys,1);
    tol = 1e-10;
    maxit = 20;
    
    [sol,er,it] = solve_with_newtons_method(f,x0,jac,tol,maxit);
    
    [sols,lm] = split_solutions(sol,fespaces_u,fespaces_p,nus,nps,indices);
    
    % plot solutions and compute errors
    
    figure(1)
    errsu = zeros(5,1);
    errsp = zeros(5,1);
    for i = 1:5
        plot_fe_fluid_function(sols{i},'U');
        hold on
        errsu(i) = compute_H1_error_velocity(fespaces_u{i},sols{i},@(x) [u1ex(x);u2ex(x)],@(x) [u1exdx(x) u1exdy(x); ...
        u2exdx(x) u2exdy(x)]);
        errsp(i) = compute_L2_error(fespaces_p{i},sols{i}.p,pex);
    end
    axis([0 1 0 1])
    hold off
    pause(0.1)
 
    errsH1u = sqrt(errsu'*errsu);
    errsL2p = sqrt(errsp'*errsp);
    
    errs = errsu + errsp;

    errH1u = [errH1u;errsH1u]
    errL2p = [errL2p;errsL2p]
    err = [err; sqrt(errs'*errs)]
end

%%
close all
n = length(h);

err = errH1u;

loglog(h(1:end),err(1:end),'.-','Markersize',10)
hold on
expo = 2;
loglog(h,h.^(expo)*min(err(1))/(h(1)^(expo)))
