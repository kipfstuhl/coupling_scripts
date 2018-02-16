clear all
close all
clc
% author: Luca Pegolotti on 10/01/2018

exsolstr = load('data/exact_solution.mat');

exsol = exsolstr.sol;

nu = 1;

fun = @(x) [0*x(1,:);0*x(2,:)];

U = 500;

dirichlet_functions = @(x) [0 0;0 0;U 0;0 0]';
neumann_functions = @(x) [0 0;0 0;0 0;0 0]';

n_elementsx = [16 23 32 45 64];
n_elementsx = [16 32 64 128];

h = 1./n_elementsx;
errH1u = [];
errL2p = [];

xline1 = 0.15;
yline1 = 0.15;

xline2 = 0.7;
yline2 = 0.35;

for nx = n_elementsx
    % create the mesh and fespaces for domain 1
    clear sol;
    xp1 = 0;
    yp1 = yline2;
    L1 = 1;
    H1 = 1-yline2;
    
    n1x = nx;
    n1y = round(n1x*H1);
    mesh1 = create_mesh(xp1,yp1,L1,H1,n1x,n1y);
    
    h_coarse = L1/n1x;
    
    bc_flags = [0 1 1 1];
    fespace1_u = create_fespace(mesh1,'P2',bc_flags);
    fespace1_p = create_fespace(mesh1,'P1',bc_flags);
    
    % create the mesh and fespaces for domain 2
    xp2 = xline2;
    yp2 = 0;
    L2 = 1-xline2;
    H2 = yline2;
    
    n2x = round(L2/h_coarse*2);
    n2y = round(H2/h_coarse*2);
    mesh2 = create_mesh(xp2,yp2,L2,H2,n2x,n2y);
    
    bc_flags = [1 1 0 0];
    fespace2_u = create_fespace(mesh2,'P2',bc_flags);
    fespace2_p = create_fespace(mesh2,'P1',bc_flags);
    
    % create the mesh and fespaces for domain 3
    xp3 = xline1;
    yp3 = 0;
    L3 = 1-xline1-L2;
    H3 = yline2;
    
    n3x = ceil(L3/h_coarse);
    n3y = ceil(H3/h_coarse);
    mesh3 = create_mesh(xp3,yp3,L3,H3,n3x,n3y);
    
    bc_flags = [1 0 0 0];
    fespace3_u = create_fespace(mesh3,'P2',bc_flags);
    fespace3_p = create_fespace(mesh3,'P1',bc_flags);
    
    % create the mesh and fespaces for domain 4
    xp4 = 0;
    yp4 = 0;
    L4 = xline1;
    H4 = yline1;
    
    n4x = round(L4/h_coarse*2);
    n4y = round(H4/h_coarse*2);
    mesh4 = create_mesh(xp4,yp4,L4,H4,n4x,n4y);
    
    bc_flags = [1 0 0 1];
    fespace4_u = create_fespace(mesh4,'P2',bc_flags);
    fespace4_p = create_fespace(mesh4,'P1',bc_flags);
    
    % create the mesh and fespaces for domain 5
    xp5 = 0;
    yp5 = yline1;
    L5 = xline1;
    H5 = 1-yline1-H1;
    
    n5x = ceil(L5/h_coarse);
    n5y = ceil(H5/h_coarse);
    mesh5 = create_mesh(xp5,yp5,L5,H5,n5x,n5y);
    
    bc_flags = [0 0 0 1];
    fespace5_u = create_fespace(mesh5,'P2',bc_flags);
    fespace5_p = create_fespace(mesh5,'P1',bc_flags);
    
    fespaces_u = {fespace1_u,fespace2_u,fespace3_u,fespace4_u,fespace5_u};
    fespaces_p = {fespace1_p,fespace2_p,fespace3_p,fespace4_p,fespace5_p};

    domain_connectivity = [1 3 3 0 3; 0 4 2 0 0; 0 0 4 2 2; 0 0 0 3 1];
    normals = [-1 0 0 0; 0 0 1 1; 0 -1 1 -1; 0 1 -1 0; 1 1 1 0];
    gausspoints = 8;
    typebasisfunctions = 'fourier';
    for bf = 1:1
        display(['bf = ', num2str(bf)])
%         nbasisfunctions = [6 5 5 4]*bf;
%        nbasisfunctions = [11 5 5 4];
        nbasisfunctions = [11 5 5 4];
        [mat,rhs,jac,nsys,nus,nps,indices] = build_coupled_system_navier_stokes(fespaces_u,fespaces_p,fun,nu,dirichlet_functions,neumann_functions,domain_connectivity,normals,nbasisfunctions,gausspoints,typebasisfunctions);
        
        % 3803
        display(['Size of the system is ', num2str(nsys)])

        % manually put a 1 on the diagonal corresponding to the degrees of
        % freedom of pressure in the bottom-left corner (corresponding to
        % domain 4)
        indices4 = indices{4};
        indexpressure = indices4(1) + 2*nus(4);
        mat = @(u) diagonalize_manually(mat(u),indexpressure);
        jac = @(u) diagonalize_manually(jac(u),indexpressure);

        % solve system with newton's method
        f = @(u) mat(u)*u-rhs;
        if (exist('sol','var'))
            updofs = 1:(2*sum(nus)+sum(nps));
            lms = sol(updofs(end)+1:end);
            lms_num = length(lms)/4;
            if (strcmp(typebasisfunctions,'polynomial'))
                lms_new = [lms(1:lms_num); zeros(2,1);
                    lms(lms_num+1:2*lms_num); zeros(2,1);
                    lms(2*lms_num+1:3*lms_num); zeros(2,1);
                    lms(3*lms_num+1:end); zeros(2,1)];
            elseif (strcmp(typebasisfunctions,'fourier'))
                lms_new = [lms(1:lms_num); zeros(4,1);
                    lms(lms_num+1:2*lms_num); zeros(4,1);
                    lms(2*lms_num+1:3*lms_num); zeros(4,1);
                    lms(3*lms_num+1:end); zeros(4,1)];
            end

            sol = sol(1:(2*sum(nus)+sum(nps)));
            x0 = [sol;lms_new];
        else
            x0 = zeros(nsys,1);
        end
        tol = 1e-8;
        maxit = 20; 

        [sol,er,it] = solve_with_newtons_method(f,x0,jac,tol,maxit);

    end
    
    [sols,lm] = split_solutions(sol,fespaces_u,fespaces_p,nus,nps,indices);
    
    intsol = interpolate_multiple_solutions(sols,exsol.fespace_u,exsol.fespace_p);
    sol1 = sols{1};
    sol2 = sols{2};
    sol4 = sols{4};
    
%     for i = 1:5
%         plot_fe_fluid_function(sols{i},'U',[0 U]);
%         hold on
%     end
    axis([0 1 0 1])
    axis square
    hold off
    pause(0.1)
    
        
    disp('=====================')
    disp(['nx = ',num2str(nx)])
    disp(['dofs vel = ',num2str(2*sum(nus))])
    disp(['dofs p = ',num2str(sum(nps))])
    disp(['lm = ',num2str(nsys-2*sum(nus)-sum(nps))])
    disp(['nsys = ',num2str(nsys)])

%     2467        9539       37507      148739
    save(['data/intsol_n',num2str(nx),'.mat'],'intsol')
    save(['data/intsoldom1_n',num2str(nx),'.mat'],'sol1')
    save(['data/intsoldom2_n',num2str(nx),'.mat'],'sol2')
    save(['data/intsoldom4_n',num2str(nx),'.mat'],'sol4')
    save(['data/nsys',num2str(nx),'.mat'],'nsys')
end

% close all
n = length(h);

err = errH1u;

loglog(h(1:end),err(1:end),'.-','Markersize',10)
hold on
expo = 2;
loglog(h,h.^(expo)*min(err(1))/(h(1)^(expo)))

function mat = diagonalize_manually(mat,row)
    mat(row,:) = 0;
    mat(row,row) = 1;
end