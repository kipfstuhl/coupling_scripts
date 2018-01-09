clear all
close all
clc

% author: Luca Pegolotti on 11/12/2017

% This script computes the error with respect to the "exact" solutiom

n_elementsx = [21 30 42 60];
h = 1./n_elementsx;
errL2u = [];
for nx = n_elementsx
    % create the mesh and fespaces for domain 1
    xp1 = 0;
    yp1 = 0.3;
    L1 = 1;
    H1 = 0.7;
    
    n1x = nx;
    n1y = floor(n1x*H1);
    mesh1 = create_mesh(xp1,yp1,L1,H1,n1x,n1y);
    
    h_coarse = L1/n1x;
    
    bc_flags = [0 1 1 1];
    fespace1_u = create_fespace(mesh1,'P2',bc_flags);
    fespace1_p = create_fespace(mesh1,'P1',bc_flags);
    
    % create the mesh and fespaces for domain 2
    xp2 = 0.7;
    yp2 = 0;
    L2 = 0.3;
    H2 = 0.3;
    
    n2x = floor(L2/h_coarse*1.5);
    n2y = floor(H2/h_coarse*1.5);
    mesh2 = create_mesh(xp2,yp2,L2,H2,n2x,n2y);
    
    bc_flags = [1 1 0 0];
    fespace2_u = create_fespace(mesh2,'P2',bc_flags);
    fespace2_p = create_fespace(mesh2,'P1',bc_flags);
    
    % create the mesh and fespaces for domain 3
    xp3 = 0.15;
    yp3 = 0;
    L3 = 0.55;
    H3 = 0.3;
    
    n3x = floor(L3/h_coarse*1.2);
    n3y = floor(H3/h_coarse*1.2);
    mesh3 = create_mesh(xp3,yp3,L3,H3,n3x,n3y);
    
    bc_flags = [1 0 0 0];
    fespace3_u = create_fespace(mesh3,'P2',bc_flags);
    fespace3_p = create_fespace(mesh3,'P1',bc_flags);
    
    % create the mesh and fespaces for domain 4
    xp4 = 0;
    yp4 = 0;
    L4 = 0.15;
    H4 = 0.15;
    
    n4x = floor(L4/h_coarse*1.5);
    n4y = floor(L4/h_coarse*1.5);
    mesh4 = create_mesh(xp4,yp2,L4,H4,n4x,n4y);
    
    bc_flags = [1 0 0 1];
    fespace4_u = create_fespace(mesh4,'P2',bc_flags);
    fespace4_p = create_fespace(mesh4,'P1',bc_flags);
    
    % create the mesh and fespaces for domain 5
    xp5 = 0;
    yp5 = 0.15;
    L5 = 0.15;
    H5 = 0.15;
    
    n5x = floor(L5/h_coarse*1.2);
    n5y = floor(H5/h_coarse*1.2);
    mesh5 = create_mesh(xp5,yp5,L5,H5,n5x,n5y);
    
    bc_flags = [0 0 0 1];
    fespace5_u = create_fespace(mesh5,'P2',bc_flags);
    fespace5_p = create_fespace(mesh5,'P1',bc_flags);
    
    fespaces_u = {fespace1_u,fespace2_u,fespace3_u,fespace4_u,fespace5_u};
    fespaces_p = {fespace1_p,fespace2_p,fespace3_p,fespace4_p,fespace5_p};
    
    % store number of degrees of freedom for one component of the velocity
    n1u = size(fespace1_u.nodes,1);
    n2u = size(fespace2_u.nodes,1);
    n3u = size(fespace3_u.nodes,1);
    n4u = size(fespace4_u.nodes,1);
    n5u = size(fespace5_u.nodes,1);
    
    % store number of degrees of freedom for the pressure
    n1p = size(fespace1_p.nodes,1);
    n2p = size(fespace2_p.nodes,1);
    n3p = size(fespace3_p.nodes,1);
    n4p = size(fespace4_p.nodes,1);
    n5p = size(fespace5_p.nodes,1);
    
    n1 = 2*n1u+n1p;
    n2 = 2*n2u+n2p;
    n3 = 2*n3u+n3p;
    n4 = 2*n4u+n4p;
    n5 = 2*n5u+n5p;
    
    indices1 = 1:n1;
    indices2 = n1+1:n1+n2;
    indices3 = n1+n2+1:n1+n2+n3;
    indices4 = n1+n2+n3+1:n1+n2+n3+n4;
    indices5 = n1+n2+n3+n4+1:n1+n2+n3+n4+n5;
    
    % define parameters and boundary conditions
    U = 1;
    f = [0;0];
    nu = 1;
    dirichlet_functions = @(x) [0 0;0 0;U*(x(2) == 1) 0;0 0]';
    neumann_functions = @(x) [0 0;0 0;0 0;0 0]';
    
    % build matrices and righ handsides for the 5 domains
    [A1,rhs1] = assembler_steady_navier_stokes(fespace1_u,fespace1_p,f,nu,dirichlet_functions,neumann_functions);
    [A2,rhs2] = assembler_steady_navier_stokes(fespace2_u,fespace2_p,f,nu,dirichlet_functions,neumann_functions);
    [A3,rhs3] = assembler_steady_navier_stokes(fespace3_u,fespace3_p,f,nu,dirichlet_functions,neumann_functions);
    [A4,rhs4] = assembler_steady_navier_stokes(fespace4_u,fespace4_p,f,nu,dirichlet_functions,neumann_functions);
    [A5,rhs5] = assembler_steady_navier_stokes(fespace5_u,fespace5_p,f,nu,dirichlet_functions,neumann_functions);
    
    
    % store number of degrees of freedom for one component the velocity
    n1u = size(fespace1_u.nodes,1);
    n2u = size(fespace2_u.nodes,1);
    n3u = size(fespace3_u.nodes,1);
    n4u = size(fespace4_u.nodes,1);
    n5u = size(fespace5_u.nodes,1);
    
    % store number of degrees of freedom for the pressure
    n1p = size(fespace1_p.nodes,1);
    n2p = size(fespace2_p.nodes,1);
    n3p = size(fespace3_p.nodes,1);
    n4p = size(fespace4_p.nodes,1);
    n5p = size(fespace5_p.nodes,1);
    
    n1 = 2*n1u+n1p;
    n2 = 2*n2u+n2p;
    n3 = 2*n3u+n3p;
    n4 = 2*n4u+n4p;
    n5 = 2*n5u+n5p;
    
    indices1 = 1:n1;
    indices2 = n1+1:n1+n2;
    indices3 = n1+n2+1:n1+n2+n3;
    indices4 = n1+n2+n3+1:n1+n2+n3+n4;
    indices5 = n1+n2+n3+n4+1:n1+n2+n3+n4+n5;
    
    A = @(u) blkdiag(A1(u(indices1)),A2(u(indices2)),A3(u(indices3)),A4(u(indices4)),A5(u(indices5)));
    
    sol_appr = zeros(n1+n2+n3+n4+n5,1);
    sollm = [];
    
    errsL2u = [];
    
    bcs_flags = [-1 0 0 0; 0 0 1 -1; 0 1 1 -1; 0 1 -1 0; 1 1 1 0]';
    
    base_freq_x = pi;
    base_freq_y = pi;
    
    B = [];
    B_t = [];
    
    for n_itx = 1:10
        for n_ity = 1:10
            gausspoints = 4;
            
            freqx = n_itx - 1;
            freqy = n_ity - 1;
            
            bf = @(x) exp(1i*x(1)*base_freq_x*freqx)*exp(1i*x(2)*base_freq_y*freqy);
            [B1,B1_t] = couple_navier_stokes_solutions_single_basis(fespaces_u,fespaces_p,bcs_flags,bf,gausspoints);
            
            B = [B;B1];
            B_t = [B_t B1_t];
        end
    end
    totalagmul = size(B,1);
    zerosp = sparse(totalagmul,totalagmul);
    mat = @(u) [A(u) B_t; B zerosp];
    
    zeroslg = zeros(totalagmul,1);
    
    % solve system with newton's method
    f = @(u) mat(u)*u-[rhs1;rhs2;rhs3;rhs4;rhs5;zeroslg];
    x0 = [sol_appr; zeros(totalagmul-size(sollm,1),1)];
    jac = @(u) [blkdiag(build_jac_navier_stokes(A1,u(indices1),fespace1_u), ...
        build_jac_navier_stokes(A2,u(indices2),fespace2_u), ...
        build_jac_navier_stokes(A3,u(indices3),fespace3_u), ...
        build_jac_navier_stokes(A4,u(indices4),fespace4_u), ...
        build_jac_navier_stokes(A5,u(indices5),fespace5_u)) B_t;
        B zerosp];
    tol = 1e-7;
    maxit = 20;
    
    [sol_apr,er,it] = solve_with_newtons_method(f,x0,jac,tol,maxit);
    
    if (it < maxit)
        sol1 = real(sol_apr(indices1));
        sol2 = real(sol_apr(indices2));
        sol3 = real(sol_apr(indices3));
        sol4 = real(sol_apr(indices4));
        sol5 = real(sol_apr(indices5));
        solm = sol_apr(n1+n2+n3+n4+n5+1:end);
        
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
        
        solstr3.u1 = sol3(1:n3u);
        solstr3.u2 = sol3(n3u+1:n3u*2);
        solstr3.p  = sol3(2*n3u+1:end);
        solstr3.fespace_u = fespace3_u;
        solstr3.fespace_p = fespace3_p;
        
        solstr4.u1 = sol4(1:n4u);
        solstr4.u2 = sol4(n4u+1:n4u*2);
        solstr4.p  = sol4(2*n4u+1:end);
        solstr4.fespace_u = fespace4_u;
        solstr4.fespace_p = fespace4_p;
        
        solstr5.u1 = sol5(1:n5u);
        solstr5.u2 = sol5(n5u+1:n5u*2);
        solstr5.p  = sol5(2*n5u+1:end);
        solstr5.fespace_u = fespace5_u;
        solstr5.fespace_p = fespace5_p;
        
        close all
        plot_fe_fluid_function(solstr1,'U',[0 U]);
        hold on
        plot_fe_fluid_function(solstr2,'U',[0 U]);
        plot_fe_fluid_function(solstr3,'U',[0 U]);
        plot_fe_fluid_function(solstr4,'U',[0 U]);
        plot_fe_fluid_function(solstr5,'U',[0 U]);
        axis([0 1 0 1])
        pause(0.01)
        
        load('data/exact_solution.mat');
        
        %     xs = 0:0.05:1;
        %     epsil = 1e-8;
        %     y1 = get_values_over_line(sol.fespace_u,sol.u2,xs,yp1,'Xpar');
        %
        %
        exact_sol_u = @(x) [evaluate_fe_function(sol.u1,sol.fespace_u,x);
            evaluate_fe_function(sol.u2,sol.fespace_u,x)];
        
        exact_sol_grad_u = @(x) [evaluate_fe_function_gradient(sol.u1,sol.fespace_u,x)';
            evaluate_fe_function_gradient(sol.u2,sol.fespace_u,x)'];
        
        ss = interpolate_multiple_solutions({solstr1,solstr2,solstr3,solstr4,solstr5},sol.fespace_u,sol.fespace_p);
        
        ss.u1 = ss.u1 - sol.u1;
        ss.u2 = ss.u2 - sol.u2
        
        plot_fe_fluid_function(ss,'U');

        pause()
        err1 = compute_H1_error_velocity(fespace1_u,solstr1,exact_sol_u,exact_sol_grad_u);
        err2 = compute_H1_error_velocity(fespace2_u,solstr2,exact_sol_u,exact_sol_grad_u);
        err3 = compute_H1_error_velocity(fespace3_u,solstr3,exact_sol_u,exact_sol_grad_u);
        err4 = compute_H1_error_velocity(fespace4_u,solstr4,exact_sol_u,exact_sol_grad_u);
        err5 = compute_H1_error_velocity(fespace5_u,solstr5,exact_sol_u,exact_sol_grad_u);
        
        errsL2u = [errsL2u; sqrt(err1^2+err2^2+err3^2+err4^2+err5^2)]
        
    else
        break;
    end
    errL2u = [errL2u;errsL2u];
end

loglog(h,errL2u)
hold on
loglog(h,h.^2*errL2u(1)/(h(1)^2));

