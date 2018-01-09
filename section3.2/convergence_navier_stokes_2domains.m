clear all
close all
clc

% author: Luca Pegolotti on 11/12/2017

u1ex = @(x) sin(x(2,:)*pi);
u2ex = @(x) exp(x(1,:))-x(1,:)*exp(1);
pex = @(x) -0.5*x(1,:)^2+0.5;

u1exdx = @(x) 0;
u1exdy = @(x) pi*cos(x(2,:)*pi);
u1exdxdx = @(x) 0;
u1exdydy = @(x) -pi^2*sin(x(2,:)*pi);

u2exdx = @(x) exp(x(1,:))-exp(1);
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

dirichlet_functions = @(x) [u1ex(x)*(x(2)==0) u2ex(x)*(x(2)==0);
                            u1ex(x)*(x(1)==1) u2ex(x)*(x(1)==1);
                            u1ex(x)*(x(2)==1) u2ex(x)*(x(2)==1);
                            u1ex(x)*(x(1)==0) u2ex(x)*(x(1)==0)]';
neumann_functions = @(x) [mu*graduex(x)*[0;-1]-pex(x)*[0;-1]*(x(2)==0), ...
                          mu*graduex(x)*[1;0]-pex(x)*[1;0]*(x(1)==1), ...
                          mu*graduex(x)*[0;1]-pex(x)*[0;1]*(x(2)==1), ...
                          mu*graduex(x)*[-1;0]-pex(x)*[-1;0]*(x(1)==0)];
n_elementsx = [8 16 32 64 128];

h = 1./n_elementsx;
errH1u = {};
errL2p = {};

for nx = n_elementsx
    % create the mesh and fespaces for domain 1
    xp1 = 0;
    yp1 = 0.5;
    L1 = 1;
    H1 = 0.5;
    
    n1x = nx;
    n1y = floor(n1x*H1);
    mesh1 = create_mesh(xp1,yp1,L1,H1,n1x,n1y);
    
    h_coarse = L1/n1x;
    
    bc_flags = [0 0 1 1];
    fespace1_u = create_fespace(mesh1,'P2',bc_flags);
    fespace1_p = create_fespace(mesh1,'P1',bc_flags);
    
    % create the mesh and fespaces for domain 2
    xp2 = 0;
    yp2 = 0;
    L2 = 1;
    H2 = 0.5;
    
    n2x = nx;
    n2y = floor(n1x*H2);
    mesh2 = create_mesh(xp2,yp2,L2,H2,n2x,n2y);
    
    h_coarse = L2/n2x;
    
    bc_flags = [1 0 0 1];
    fespace2_u = create_fespace(mesh2,'P2',bc_flags);
    fespace2_p = create_fespace(mesh2,'P1',bc_flags);    
    
    % store number of degrees of freedom for one component of the velocity
    n1u = size(fespace1_u.nodes,1);
    n2u = size(fespace2_u.nodes,1);
    
    % store number of degrees of freedom for the pressure
    n1p = size(fespace1_p.nodes,1);
    n2p = size(fespace2_p.nodes,1);

    n1 = 2*n1u+n1p;
    n2 = 2*n2u+n2p;
    
    indices1 = 1:n1;
    indices2 = n1+1:n1+n2;

    % build matrices and righ handsides for the 5 domains
    [A1,rhs1] = assembler_steady_navier_stokes(fespace1_u,fespace1_p,fun,nu,dirichlet_functions,neumann_functions);
    [A2,rhs2] = assembler_steady_navier_stokes(fespace2_u,fespace2_p,fun,nu,dirichlet_functions,neumann_functions);
    
    % store number of degrees of freedom for one component the velocity
    n1u = size(fespace1_u.nodes,1);
    n2u = size(fespace2_u.nodes,1);

    % store number of degrees of freedom for the pressure
    n1p = size(fespace1_p.nodes,1);
    n2p = size(fespace2_p.nodes,1);
    
    n1 = 2*n1u+n1p;
    n2 = 2*n2u+n2p;
 
    indices1 = 1:n1;
    indices2 = n1+1:n1+n2;

    A = @(u) blkdiag(A1(u(indices1)),A2(u(indices2)));
    
    sol_appr = zeros(n1+n2,1);
    sollm = [];
    
    errsH1u = [];
    errsL2p = [];
    for n_it = 1:10
        
        n_iterations = n_it;
        gausspoints = 4;
        
        % interface 1
        fespaces_u = {fespace1_u, fespace2_u};
        fespaces_p = {fespace1_p, fespace2_p};
        bcs_flags = [1 0 0 0; 0 0 1 0]';
        base_freq = 1/L1;
        label = 'xpar';
        n_iterations1 = n_iterations;
        [blocks_interface1,blocks_interface1_t] = couple_navier_stokes_solutions(fespaces_u,fespaces_p,bcs_flags,base_freq,n_iterations1,label,gausspoints);
        nlagmul1 = (n_iterations1-1)*4 + 2;
        
        
        B1 = [-blocks_interface1{1} blocks_interface1{2}];

        B = [B1];
        
        B1_t = [-blocks_interface1_t{1};blocks_interface1_t{2}];

        B_t = [B1_t];
        
        totalagmul = nlagmul1;
        
        zerosp = sparse(totalagmul,totalagmul);
        mat = @(u) [A(u) B_t; B zerosp];
        
        zeroslg = zeros(totalagmul,1);
        
        % solve system with newton's method
        f = @(u) mat(u)*u-[rhs1;rhs2;zeroslg];
        x0 = [sol_appr; zeros(totalagmul-size(sollm,1),1)];
        jac = @(u) [blkdiag(build_jac_navier_stokes(A1,u(indices1),fespace1_u), ...
            build_jac_navier_stokes(A2,u(indices2),fespace2_u)) B_t;
            B zerosp];
        tol = 1e-7;
        maxit = 20;
        
        [sol_apr,er,it] = solve_with_newtons_method(f,x0,jac,tol,maxit);
        
        if (it < maxit)
            sol1 = sol_apr(indices1);
            sol2 = sol_apr(indices2);
          
            solm = sol_apr(n1+n2+1:end);
            
            solstr1.u1 = sol1(1:n1u);
            solstr1.u2 = sol1(n1u+1:n1u*2);
            solstr1.p  = sol1(2*n1u+1:end);
            solstr1.fespace_u = fespace1_u;
            solstr1.fespace_p = fespace1_p;
            
            solstr2.u1 = sol2(1:n2u);
            solstr2.u2 = sol2(n2u+1:n2u*2);
            solstr2.p  = sol2(2*n2u+1:end);
            solstr2.fespace_u = fespace2_u;
            solstr2.fespace_p = fespace2_p;
                
            close all
            plot_fe_fluid_function(solstr1,'U',[0 1]);
            hold on
            plot_fe_fluid_function(solstr2,'U',[0 1]);

            axis([0 1 0 1])
            pause(0.01)
            
            err1 = compute_H1_error_velocity(fespace1_u,solstr1,@(x) [u1ex(x);u2ex(x)],@(x) [u1exdx(x) u1exdy(x); ...
        u2exdx(x) u2exdy(x)]);
            err2 = compute_H1_error_velocity(fespace2_u,solstr2,@(x) [u1ex(x);u2ex(x)],@(x) [u1exdx(x) u1exdy(x); ...
        u2exdx(x) u2exdy(x)]);
    
            err1p = compute_L2_error(fespace1_p,solstr1.p,pex);
            err2p = compute_L2_error(fespace2_p,solstr2.p,pex);

            
            errsH1u = [errsH1u; sqrt(err1^2+err2^2)]
            errsL2p = [errsL2p; sqrt(err1p^2+err2p^2)]

        else
            break;
        end
    end
    errH1u{end+1} = errsH1u;
    errL2p{end+1} = errsL2p;
end

% loglog(h,errL2u)
% hold on
% loglog(h,h.^3*errL2u(1)/(h(1)^3));
%%
n = length(h);

err = [];

for i = 1:n
    err = [err;min(errL2p{i}+errH1u{i})];
end

loglog(h,err)
hold on
loglog(h,h.^2*min(err(1))/(h(1)^1))
