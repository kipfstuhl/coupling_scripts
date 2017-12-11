clear all
close all
clc

% author: Luca Pegolotti on 1/12/2017

% This script computes the convergence order obtained on non conforming meshes
% presented in Section 3. The Poisson problem is solved on two subdomains:
% (0,0.5)x(0,1) and (0.5,0)x(0,1). The basis functions on the interface are
% Fourier basis functions.

% NOTE: the script is actually the same as
% compute_order_convergence_P2P2_confmesh.m exept for little details. For
% this reason, the code is cleared from unnecessary comments.

set(0,'defaulttextinterpreter','latex')

run load_exact_solution_and_f.m

dir_functions = @(x) [0;0;0;0];
neu_functions = @(x) [0;0;0;0];

N = [20 28 40 56 80 114 160];

typerror = 'H1';

betah = [];

% we need high order quadrature for the non-conforming mesh (we invite the
% user to try with e.g. quadrature_order = 3 or quadrature_order = 5)
quadrature_order = 7;

for n_elements = N    
    n1x = n_elements/2;
    n2x = n_elements/2;
    n1y = n_elements;
    n2y = n_elements;
    
    xp1 = 0;
    yp1 = 0;
    L1 = 0.5;
    H1 = 1;
    
    mesh1 = create_mesh(xp1,yp1,L1,H1,n1x,n1y);
    
    xp2 = L1;
    yp2 = 0;
    L2 = 0.5;
    H2 = 1;
    
    mesh2 = create_mesh(xp2,yp2,L2,H2,n2x,n2y);
    
    bc1 = [1 0 1 1];
    fespace1 = create_fespace(mesh1,'P2',bc1);
    
    [A1,rhs1] = assembler_poisson(fespace1,fun,mu,dir_functions,neu_functions);
    
    % create also mass matrix for first problem (needed to compute norm
    % matrix)
    M1 = assemble_mass(fespace1);
    
    bc2 = [1 1 1 0];
    fespace2 = create_fespace(mesh2,'P2',bc2);
    
    [A2,rhs2] = assembler_poisson(fespace2,fun,mu,dir_functions,neu_functions);
    M2 = assemble_mass(fespace2);
    
    n1 = size(A1,1);
    n2 = size(A2,1);
    
    indices1 = 1:n1;
    indices2 = n1+1:n1+n2;
    
    n_iterations = 16;
    
    B1 = [];
    B2 = [];
    
    betahs = [];
    Mlambda = [];
    
    for i = 1:n_iterations
        disp(['n_elements = ', num2str(n_elements) ...
              ', iteration n = ', num2str(i)]);
        
        b1 = zeros(n1,1);
        b2 = zeros(n2,1);
        
        freq = i - 1;

        if (i == 1)
            b1 = apply_neumann_bc(fespace1,b1,@(x) [0;1;0;0]);
            
            B1 = [B1;b1'];
            
            b2 = apply_neumann_bc(fespace2,b2,@(x) [0;0;0;1]);
            
            B2 = [B2;b2'];
            Mlambda = 1;
        else
            b1 = apply_neumann_bc(fespace1,b1,@(x) [0;sin(x(2) * pi * freq);0;0],quadrature_order);
            B1 = [B1;b1'];
            
            b1 = b1*0;
            
            b1 = apply_neumann_bc(fespace1,b1,@(x) [0;cos(x(2) * pi * freq);0;0],quadrature_order);
            B1 = [B1;b1'];
            
            b2 = apply_neumann_bc(fespace2,b2,@(x) [0;0;0;sin(x(2) * pi * freq)],quadrature_order);
            B2 = [B2;b2'];
            
            b2 = b2*0;
            
            b2 = apply_neumann_bc(fespace2,b2,@(x) [0;0;0;cos(x(2) * pi * freq)],quadrature_order);
            B2 = [B2;b2'];
            
            Mlambda = [Mlambda zeros(2*(i-2)+1,2);
                       zeros(2,2*(i-2)+1) zeros(2,2)];
                   
            % fill mass matrix for Lambda space
            Mlambda(end-1,1) = integral(@(x) sin(x*pi*freq),0,1);
            Mlambda(end,1)   = integral(@(x) cos(x*pi*freq),0,1);
            for j = 2:2:2*(i-1)
                freq2 = j/2;
                Mlambda(j,end-1) = integral(@(x) sin(x*pi*freq).*sin(x*pi*freq2),0,1);
                Mlambda(end-1,j) = Mlambda(j,end-1);
            end
            
            for j = 3:2:2*(i-1)
                freq2 = (j-1)/2;
                Mlambda(j,end-1) = integral(@(x) sin(2*x*pi*freq).*cos(2*x*pi*freq2),0,1);
                Mlambda(end-1,j) = Mlambda(j,end-1);
            end
            
            for j = 2:2:2*(i-1)+1
                freq2 = j/2;
                Mlambda(j,end) = integral(@(x) cos(x*pi*freq).*sin(x*pi*freq2),0,1);
                Mlambda(end,j) = Mlambda(j,end);
            end
            
            for j = 3:2:2*(i-1)+1
                freq2 = (j-1)/2;
                Mlambda(j,end) = integral(@(x) cos(x*pi*freq).*cos(x*pi*freq2),0,1);
                Mlambda(end,j) = Mlambda(j,end);
            end
        end        
        n3 = size(B1,1);
        
        leftmat = [A1+M1 sparse(n1,n2) -B1'; ...
                   sparse(n2,n1) A2+M2  B2'; ...
                   -B1 B2 sparse(n3,n3)];
        
        rightmat = sparse(size(leftmat,1),size(leftmat,1));
        rightmat(end-n3+1:end,end-n3+1:end) = Mlambda;
        
        e = -eigs(leftmat,rightmat,1,'SM','tolerance',1e-16);
        betahs = [betahs;sqrt(e)]
    end
    betah = [betah betahs];
end

if (exist('data','dir') == 0)
    disp('error')
    mkdir data;
end
save('data/betah.mat','betah');
