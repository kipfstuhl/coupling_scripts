clear all
close all
clc

% author: Luca Pegolotti on 11/12/2017
run load_exact_solution_and_f.m
dirichlet_functions = @(x) [0;0;0;0];
neumann_functions = @(x) [x(1,:)*0;uexdx(x(1,:),x(2,:)).*(x(1,:)==1);x(1,:)*0;x(1,:)*0];

n_elementsx = [16 32 64 128];

h = 1./n_elementsx;
errH1u = [];
errL2p = [];
err = [];

xline = 0.7;
yline = 0.5;
for nx = n_elementsx
    % create the mesh and fespaces for domain 1
    xp1 = 0;
    yp1 = yline;
    L1 = 1;
    H1 = 1-yline;
    
    n1x = nx;
    n1y = floor(n1x*H1);
    mesh1 = create_mesh(xp1,yp1,L1,H1,n1x,n1y);
    
    h_coarse = L1/n1x;
    
    bc_flags = [0 1 1 1];
    fespace1 = create_fespace(mesh1,'P2',bc_flags);
    
    % create the mesh and fespaces for domain 2
    xp2 = xline;
    yp2 = 0;
    L2 = 1-xline;
    H2 = yline;
    
    n2x = floor(nx*L2);
    n2y = floor(nx*H2);
    mesh2 = create_mesh(xp2,yp2,L2,H2,n2x,n2y);
    
    bc_flags = [1 1 0 0];
    fespace2 = create_fespace(mesh2,'P2',bc_flags);
    
    % create the mesh and fespaces for domain 3
    xp3 = 0;
    yp3 = 0;
    L3 = xline;
    H3 = yline;
    
    n3x = floor(nx*L3);
    n3y = floor(nx*H3);
    mesh3 = create_mesh(xp3,yp3,L3,H3,n3x,n3y);
    
    bc_flags = [1 0 0 1];
    fespace3 = create_fespace(mesh3,'P2',bc_flags);
    
    
    fespaces = {fespace1,fespace2,fespace3};
    
    figure(100)
    multimesh = {mesh1,mesh2,mesh3};
    draw_multimesh(multimesh);
    
    pause()
    
    % store number of degrees of freedom for one component of the velocity
    n1 = size(fespace1.nodes,1);
    n2 = size(fespace2.nodes,1);
    n3 = size(fespace3.nodes,1);
    
    indices1 = 1:n1;
    indices2 = n1+1:n1+n2;
    indices3 = n1+n2+1:n1+n2+n3;
    
    % build matrices and righ handsides for the 5 domains
    [A1,rhs1] = assembler_poisson(fespace1,fun,mu,dirichlet_functions,neumann_functions);
    [A2,rhs2] = assembler_poisson(fespace2,fun,mu,dirichlet_functions,neumann_functions);
    [A3,rhs3] = assembler_poisson(fespace3,fun,mu,dirichlet_functions,neumann_functions);
    
    A = blkdiag(A1,A2,A3);
    
    sol_appr = zeros(n1+n2+n3,1);
    sollm = [];
    
    errsL2u = [];
    
    bcs_flags = [-1 0 0 0; 0 0 1 -1; 0 1 1 0]';
    
    
    B = [];
    B_t = [];
    
    gausspoints = 40;
    
%     bf = @(x) (x(2) == 0.5);
%     [B1,B1_t] = couple_poisson_solutions(fespaces,bcs_flags,bf,gausspoints);
%     B = [B;B1];
%     B_t = [B_t B1_t];
%     
%     bf = @(x) (x(1) == 0.5);
%     [B1,B1_t] = couple_poisson_solutions(fespaces,bcs_flags,bf,gausspoints);
%     B = [B;B1];
%     B_t = [B_t B1_t];
%     
    for n_itx = 0:7
        bf = @(x) (abs(x(2) -yline)<1e-6)*x(1).^n_itx;
        [B1,B1_t] = couple_poisson_solutions(fespaces,bcs_flags,bf,gausspoints);
        B = [B;B1];
        B_t = [B_t B1_t];
    end
    
    for n_ity = 0:7
        bf = @(x) (abs(x(1)-xline)<1e-6)*x(2).^n_ity;
        [B1,B1_t] = couple_poisson_solutions(fespaces,bcs_flags,bf,gausspoints);
        B = [B;B1];
        B_t = [B_t B1_t];
    end
    totalagmul = size(B,1);
    zerosp = sparse(totalagmul,totalagmul);
    mat = [A B_t; B zerosp];
    
    %sum(B_t,1)
    
    zeroslg = zeros(totalagmul,1);
    sol_apr = mat\[rhs1;rhs2;rhs3;zeroslg];
    
    sol1 = real(sol_apr(indices1));
    sol2 = real(sol_apr(indices2));
    sol3 = real(sol_apr(indices3));
    
    solm = sol_apr(n1+n2+n3+1:end);
    figure(1)
    plot_fe_function(sol1,fespace1);
    hold on
    plot_fe_function(sol2,fespace2);
    plot_fe_function(sol3,fespace3);
    hold off
    
    pause(0.1)
    
    % plot the derivatives at the interface and the lagrange multiplier
    figure(2)
    box
    hold off
    % we compute the derivatives at x = 0.5
    y0 = yline;
    
    % epsilon used for the numerical approximation via finite differences
    % of the derivatives
    epsil = 1e-9;
    
    % approximate derivative on the right subdomain
    [xx1,u1] = get_values_over_line(fespace1,sol1,100,y0+epsil,'Xpar');
    [~ ,u2] = get_values_over_line(fespace1,sol1,100,y0,'Xpar');
    hold on
    
    plot(xx1,(u2-u1)/epsil,'r')
    
    % approximate derivative on the right subdomain
    [xx2,u1] = get_values_over_line(fespace2,sol2,100,y0,'Xpar');
    [~ ,u2] = get_values_over_line(fespace2,sol2,100,y0-epsil,'Xpar');
    
    plot(xx2,(u2-u1)/epsil,'b')
    
    % approximate derivative on the right subdomain
    [xx3,u1] = get_values_over_line(fespace3,sol3,100,y0,'Xpar');
    [~ ,u2] = get_values_over_line(fespace3,sol3,100,y0-epsil,'Xpar');
    
    plot(xx3,(u2-u1)/epsil,'y')
    
    % plot the exact derivative
    plot(xx1,-uexdy(xx1,y0*xx1.^0),'--k')
    
    % plot the Lagrange multiplier
    coeffs = sol_apr(n1+n2+n3+1:end);
    urec = xx1*0 - coeffs(1);
    
    coeffs = coeffs(2:end);
    
    for j = 1:n_itx
        urec = urec + xx1.^(j)*coeffs(j);
    end
    
    plot(xx1,urec,'color',[0 147 23]/255)
    hold off
    
    % plot the derivatives at the interface and the lagrange multiplier
    figure(3)
    box
    hold off
    % we compute the derivatives at x = 0.5
    x0 = xline;
    
    % epsilon used for the numerical approximation via finite differences
    % of the derivatives
    epsil = 1e-9;
    
    % approximate derivative on the right subdomain
    [x2,u1] = get_values_over_line(fespace2,sol2,100,x0+epsil,'Ypar');
    [~ ,u2] = get_values_over_line(fespace2,sol2,100,x0,'Ypar');
    hold on
    
    plot(x2,(u2-u1)/epsil,'b')
    
    % approximate derivative on the right subdomain
    [x3,u1] = get_values_over_line(fespace3,sol3,100,x0,'Ypar');
    [~ ,u2] = get_values_over_line(fespace3,sol3,100,x0-epsil,'Ypar');
       
    plot(x3,(u2-u1)/epsil,'r')
    
    % plot the exact derivative
    plot(x2,-uexdx(x0*x2.^0,x2),'--k')
    
    % plot the Lagrange multiplier
    coeffs = sol_apr(n1+n2+n3+1:end);
    urec = x2*0 - coeffs(2+n_itx);
    
    coeffs = coeffs(n_itx+3:end);
    
    for j = 1:n_ity
        urec = urec + x2.^(j)*coeffs(j);
    end
    
    plot(x2,urec,'color',[0 147 23]/255)
    hold off
    
    pause()
    close all
    
    err1 = compute_L2_error(fespace1,sol1,uex);
    err2 = compute_L2_error(fespace2,sol2,uex);
    err3 = compute_L2_error(fespace3,sol3,uex);
    
    
    err = [err; sqrt(err1^2 + err2^2 + err3^2)]
end

% loglog(h,errL2u)
% hold on
% loglog(h,h.^3*errL2u(1)/(h(1)^3));
%%
n = length(h);

%err = errH1u + errL2p;

loglog(h,err)
hold on
expo = 3;
loglog(h,h.^(expo)*min(err(1))/(h(1)^(expo)))
