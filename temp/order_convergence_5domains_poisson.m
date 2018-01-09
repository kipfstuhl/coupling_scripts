clear all
close all
clc

% author: Luca Pegolotti on 11/12/2017
run load_exact_solution_and_f.m
dirichlet_functions = @(x) [0;0;0;0];
neumann_functions = @(x) [0;0;0;0];

n_elementsx = [16 32 64 128];

h = 1./n_elementsx;
errH1u = [];
errL2p = [];
err = [];
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
    
    bc_flags = [0 1 1 1];
    fespace1 = create_fespace(mesh1,'P2',bc_flags);
    
    % create the mesh and fespaces for domain 2
    xp2 = 0.75;
    yp2 = 0;
    L2 = 0.25;
    H2 = 0.5;
    
    n2x = nx/4;
    n2y = nx/2;
    mesh2 = create_mesh(xp2,yp2,L2,H2,n2x,n2y);
    
    bc_flags = [1 1 0 0];
    fespace2 = create_fespace(mesh2,'P2',bc_flags);
    
    % create the mesh and fespaces for domain 3
    xp3 = 0.25;
    yp3 = 0;
    L3 = 0.5;
    H3 = 0.5;
    
    n3x = nx/2;
    n3y = nx/2;
    mesh3 = create_mesh(xp3,yp3,L3,H3,n3x,n3y);
    
    bc_flags = [1 0 0 0];
    fespace3 = create_fespace(mesh3,'P2',bc_flags);
    
    % create the mesh and fespaces for domain 4
    xp4 = 0;
    yp4 = 0;
    L4 = 0.25;
    H4 = 0.25;
    
    n4x = nx/4;
    n4y = nx/4;
    mesh4 = create_mesh(xp4,yp4,L4,H4,n4x,n4y);
    
    bc_flags = [1 0 0 1];
    fespace4 = create_fespace(mesh4,'P2',bc_flags);
    
    % create the mesh and fespaces for domain 5
    xp5 = 0;
    yp5 = 0.25;
    L5 = 0.25;
    H5 = 0.25;
    
    n5x = nx/4;
    n5y = nx/4;
    mesh5 = create_mesh(xp5,yp5,L5,H5,n5x,n5y);
    
    bc_flags = [0 0 0 1];
    fespace5 = create_fespace(mesh5,'P2',bc_flags);
    
    fespaces = {fespace1,fespace2,fespace3,fespace4,fespace5};
    
    %     meshes = {};
    %     meshes{end+1} = mesh1;
    %     meshes{end+1} = mesh2;
    %     meshes{end+1} = mesh3;
    %     meshes{end+1} = mesh4;
    %     meshes{end+1} = mesh5;
    %
    %     draw_multimesh(meshes)
    %     return;
    % store number of degrees of freedom for one component of the velocity
    n1 = size(fespace1.nodes,1);
    n2 = size(fespace2.nodes,1);
    n3 = size(fespace3.nodes,1);
    n4 = size(fespace4.nodes,1);
    n5 = size(fespace5.nodes,1);
    
    
    indices1 = 1:n1;
    indices2 = n1+1:n1+n2;
    indices3 = n1+n2+1:n1+n2+n3;
    indices4 = n1+n2+n3+1:n1+n2+n3+n4;
    indices5 = n1+n2+n3+n4+1:n1+n2+n3+n4+n5;
    
    % build matrices and righ handsides for the 5 domains
    [A1,rhs1] = assembler_poisson(fespace1,fun,mu,dirichlet_functions,neumann_functions);
    [A2,rhs2] = assembler_poisson(fespace2,fun,mu,dirichlet_functions,neumann_functions);
    [A3,rhs3] = assembler_poisson(fespace3,fun,mu,dirichlet_functions,neumann_functions);
    [A4,rhs4] = assembler_poisson(fespace4,fun,mu,dirichlet_functions,neumann_functions);
    [A5,rhs5] = assembler_poisson(fespace5,fun,mu,dirichlet_functions,neumann_functions);
    
    
    A = blkdiag(A1,A2,A3,A4,A5);
    
    sol_appr = zeros(n1+n2+n3+n4+n5,1);
    sollm = [];
    
    errsL2u = [];
    
    bcs_flags = [-1 0 0 0; 0 0 1 -1; 0 1 1 -1; 0 1 -1 0; 1 1 1 0]';
    
    base_freq_x = pi;
    base_freq_y = pi;
    
    multx1 = 0.47;
    multx2 = 0.33;
    
    multy1 = 0.47;
    multy2 = 0.33;
    
    B = [];
    B_t = [];
    
    %     for n_itx = 1:6
    %         for n_ity = 1:6
    %             gausspoints = 4;
    %
    %             freqx = n_itx - 1;
    %             freqy = n_ity - 1;
    %
    %
    %             bf = @(x) exp(1i*x(1)*base_freq_x*freqx)*exp(1i*x(2)*base_freq_y*freqy);
    %             [B1,B1_t] = couple_navier_stokes_solutions_single_basis(fespaces_u,fespaces_p,bcs_flags,bf,gausspoints);
    %
    %             B = [B;B1];
    %             B_t = [B_t B1_t];
    %         end
    %     end
    gausspoints = 4;
    
    for n_itx = 1:20
        freqx = n_itx - 1;
        
        if (freqx ~= 0)
            bf = @(x) sin(x(1)*base_freq_x*freqx)*sin(x(2)*pi*multy1);
            [B1,B1_t] = couple_poisson_solutions(fespaces,bcs_flags,bf,gausspoints);
            B = [B;B1];
            B_t = [B_t B1_t];
            
            bf = @(x) cos(x(1)*base_freq_x*freqx)*sin(x(2)*pi*multy1);
            [B1,B1_t] = couple_poisson_solutions(fespaces,bcs_flags,bf,gausspoints);
            B = [B;B1];
            B_t = [B_t B1_t];
            
            bf = @(x) sin(x(1)*base_freq_x*freqx)*sin(x(2)*pi*multy2);
            [B1,B1_t] = couple_poisson_solutions(fespaces,bcs_flags,bf,gausspoints);
            B = [B;B1];
            B_t = [B_t B1_t];
            
            bf = @(x) cos(x(1)*base_freq_x*freqx)*sin(x(2)*pi*multy2);
            [B1,B1_t] = couple_poisson_solutions(fespaces,bcs_flags,bf,gausspoints);
            B = [B;B1];
            B_t = [B_t B1_t];
        end
    end
    
    for n_ity = 1:20
        freqy = n_ity - 1;
        
        if (freqy ~= 0)
            bf = @(x) sin(x(2)*base_freq_y*freqy)*sin(x(1)*pi*multx1);
            [B1,B1_t] = couple_poisson_solutions(fespaces,bcs_flags,bf,gausspoints);
            B = [B;B1];
            B_t = [B_t B1_t];
            
            bf = @(x) cos(x(2)*base_freq_y*freqy)*sin(x(1)*pi*multx1);
            [B1,B1_t] = couple_poisson_solutions(fespaces,bcs_flags,bf,gausspoints);
            B = [B;B1];
            B_t = [B_t B1_t];
            
            bf = @(x) sin(x(2)*base_freq_y*freqy)*sin(x(1)*pi*multx2);
            [B1,B1_t] = couple_poisson_solutions(fespaces,bcs_flags,bf,gausspoints);
            B = [B;B1];
            B_t = [B_t B1_t];
            
            bf = @(x) cos(x(2)*base_freq_y*freqy)*sin(x(1)*pi*multx2);
            [B1,B1_t] = couple_poisson_solutions(fespaces,bcs_flags,bf,gausspoints);
            B = [B;B1];
            B_t = [B_t B1_t];
        else
            bf = @(x) 1;
            [B1,B1_t] = couple_poisson_solutions(fespaces,bcs_flags,bf,gausspoints);
            B = [B;B1];
            B_t = [B_t B1_t];
        end
    end
    totalagmul = size(B,1);
    zerosp = sparse(totalagmul,totalagmul);
    mat = [A B_t; B zerosp];
    
    zeroslg = zeros(totalagmul,1);
    sol_apr = mat\[rhs1;rhs2;rhs3;rhs4;rhs5;zeroslg];
    
    sol1 = real(sol_apr(indices1));
    sol2 = real(sol_apr(indices2));
    sol3 = real(sol_apr(indices3));
    sol4 = real(sol_apr(indices4));
    sol5 = real(sol_apr(indices5));
    solm = sol_apr(n1+n2+n3+n4+n5+1:end);
    
    plot_fe_function(sol1,fespace1);
    hold on
    plot_fe_function(sol2,fespace2);
    plot_fe_function(sol3,fespace3);
    plot_fe_function(sol4,fespace4);
    plot_fe_function(sol5,fespace5);
    hold off
    
    pause(0.1)
    
    
    err1 = compute_L2_error(fespace1,sol1,uex);
    err2 = compute_L2_error(fespace2,sol2,uex);
    err3 = compute_L2_error(fespace3,sol3,uex);
    err4 = compute_L2_error(fespace4,sol4,uex);
    err5 = compute_L2_error(fespace5,sol5,uex);
    
    
    err = [err; sqrt(err1^2 + err2^2 + err3^2 + err4^2 + err5^2)]
end

% loglog(h,errL2u)
% hold on
% loglog(h,h.^3*errL2u(1)/(h(1)^3));
%%
n = length(h);

%err = errH1u + errL2p;

loglog(h,err)
hold on
expo = 1.4;
loglog(h,h.^(expo)*min(err(1))/(h(1)^(expo)))
