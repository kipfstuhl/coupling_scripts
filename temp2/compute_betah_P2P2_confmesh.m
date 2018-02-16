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
N = [20 28 40 56 80];

colors = [1 0 0;
         0 1 0;
         0 0 1;
         0.5 0.5 0;
         0 0.5 0.5];
%N = [20 28];
typerror = 'H1';

betah_l2l2 = {};
betah_h1l2 = {};
betah_h1h1 = {};

nbfsh = {};
nith = {};
condh = {};

% we need high order quadrature for the non-conforming mesh (we invite the
% user to try with e.g. quadrature_order = 3 or quadrature_order = 5)
quadrature_order = 3;
count = 0;
for n_elements = N
    count = count + 1;
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
    %M1 = apply_dirichlet_bc_matrix(assemble_mass(fespace1),fespace1,1);
    D1 = assemble_stiffness(1,fespace1);
    M1 = assemble_mass(fespace1);
    
    bc2 = [1 1 1 0];
    fespace2 = create_fespace(mesh2,'P2',bc2);
    
    [A2,rhs2] = assembler_poisson(fespace2,fun,mu,dir_functions,neu_functions);
    %M2 = apply_dirichlet_bc_matrix(assemble_mass(fespace2),fespace2,1);
    D2 = assemble_stiffness(1,fespace2);
    M2 = assemble_mass(fespace2);
    
    n1 = size(A1,1);
    n2 = size(A2,1);
    
    indices1 = 1:n1;
    indices2 = n1+1:n1+n2;
    
    n_iterations = 400;
    
    B1 = [];
    B2 = [];
    
    B1_t = [];
    B2_t = [];
    
    betahs_l2l2 = [];
    betahs_h1l2 = [];
    betahs_h1h1 = [];
    nbfs = [];
    Mlambda = [];
    bfs = {};
    bfsd = {};
    conds = [];
    
    
    for i = 1:n_iterations
        disp(['n_elements = ', num2str(n_elements) ...
            ', iteration n = ', num2str(i)]);
        
        b1 = zeros(n1,1);
        b2 = zeros(n2,1);
        
        freq = (i - 1);
        
        if (i == 1)
            b1 = apply_neumann_bc(b1,fespace1,@(x) [0;1;0;0]);
            
            B1 = [B1;b1'];
            
            B1_t = [B1_t apply_dirichlet_bc_rhs(b1,fespace1,@(x) [0;0;0;0])];
            
            b2 = apply_neumann_bc(b2,fespace2,@(x) [0;0;0;1]);
            
            B2 = [B2;b2'];
            B2_t = [B2_t apply_dirichlet_bc_rhs(b2,fespace2,@(x) [0;0;0;0])];
            
            bfs{end+1} = @(x) x.^0;
            bfsd{end+1} = @(x) 0*x.^0;
        else
            b1 = apply_neumann_bc(b1,fespace1,@(x) [0;sin(x(2) * pi * freq);0;0],quadrature_order);
            B1 = [B1;b1'];
            
            B1_t = [B1_t apply_dirichlet_bc_rhs(b1,fespace1,@(x) [0;0;0;0])];
            
            
            b1 = b1*0;
            
            bfs{end+1} = @(x) sin(x * pi * freq);
            bfsd{end+1} = @(x) pi * freq*cos(x * pi * freq);
            
            b1 = apply_neumann_bc(b1,fespace1,@(x) [0;cos(x(2) * pi * freq);0;0],quadrature_order);
            B1 = [B1;b1'];
            
            B1_t = [B1_t apply_dirichlet_bc_rhs(b1,fespace1,@(x) [0;0;0;0])];
            
            
            bfs{end+1} = @(x) cos(x * pi * freq);
            bfsd{end+1} = @(x) -pi * freq*sin(x * pi * freq);
            
            b2 = apply_neumann_bc(b2,fespace2,@(x) [0;0;0;sin(x(2) * pi * freq)],quadrature_order);
            B2 = [B2;b2'];
            
            B2_t = [B2_t apply_dirichlet_bc_rhs(b2,fespace2,@(x) [0;0;0;0])];
            
            b2 = b2*0;
            
            b2 = apply_neumann_bc(b2,fespace2,@(x) [0;0;0;cos(x(2) * pi * freq)],quadrature_order);
            B2 = [B2;b2'];
            
            B2_t = [B2_t apply_dirichlet_bc_rhs(b2,fespace2,@(x) [0;0;0;0])];
        end
        
        %try
            if (i > 0 || floor(log2(i)) == log2(i))
                lambdas = size(bfs,2);
                Mlambda = zeros(lambdas);
                for k1 = 1:lambdas
                    for k2 = 1:lambdas
                        Mlambda(k1,k2) = integral(@(x) bfs{k1}(x).*bfs{k2}(x),0,1);
                    end
                end
                n3 = size(B1,1);
                
                leftmat = [M1 sparse(n1,n2) -B1'; ...
                    sparse(n2,n1) M2  B2'; ...
                    -B1 B2 sparse(n3,n3)];
                
 
                B = [-B1 B2];
                % e = abs(invitr(leftmat\rightmat,1000,1e-15));
                mat = B*([M1 sparse(n1,n2);sparse(n2,n1) M2]\B');
                if (size(mat,1) == 1)
                    e = mat/Mlambda;
                else
                    e = eigs(B*([M1 sparse(n1,n2);sparse(n2,n1) M2]\B'),Mlambda,1,'SM')
                end
                % e = abs(eigs(leftmat,rightmat,1,'SM','tolerance',1e-30));
                
                figure(1)
                title('L2-L2')
                semilogy(i,sqrt(e),'.','color',colors(count,:),'Markersize',10)
                hold on
%                 betahs_l2l2 = [betahs_l2l2;sqrt(e)]
%                 nbfs = [nbfs;size(bfs,2)];
%                 mat = [A1 sparse(n1,n2) -B1_t; ...
%                     sparse(n2,n1) A2 B2_t; ...
%                     -B1 B2 sparse(n3,n3)];
%                 conds = [conds;condest(mat)]
                
%                 B = [-B1 B2];
%                 Bt = sparse([-B1_t;B2_t]);
%                 % e = abs(invitr(leftmat\rightmat,1000,1e-15));
%                 mmm = B*([M1 sparse(n1,n2);sparse(n2,n1) M2]\B');
%                 e = abs(eigs(B*([M1 sparse(n1,n2);sparse(n2,n1) M2]\B'),Mlambda,1,'SM'));
%                 pause()
%                 %e = abs(eigs(B*([M1 sparse(n1,n2);sparse(n2,n1) M2]\B'),Mlambda,1,'SM','tolerance',1e-30));
%                 % e = abs(eigs(leftmat,rightmat,1,'SM','tolerance',1e-30));
%                 %e = invitr(B*([M1 sparse(n1,n2);sparse(n2,n1) M2]\B'),Mlambda,1000,1e-16);
%                 figure(1)
%                 title('H1-L2')
%                 semilogy(i,sqrt(e),'.','color',colors(count,:),'Markersize',10)
%                 hold on
%                 betahs_h1l2 = [betahs_h1l2;sqrt(e)]
%                 if (betahs_h1l2(end) < 1e-9)
%                     break
%                 end
                

%                 
%                 for k1 = 1:lambdas
%                     for k2 = 1:lambdas
%                         Mlambda(k1,k2) = Mlambda(k1,k2)+ integral(@(x) bfsd{k1}(x).*bfsd{k2}(x),0,1);
%                     end
%                 end
%                 
%                 leftmat = [M1+D1 sparse(n1,n2) -B1'; ...
%                     sparse(n2,n1) M2+D2  B2'; ...
%                     -B1 B2 sparse(n3,n3)];
%                 
%                 rightmat(end-n3+1:end,end-n3+1:end) = Mlambda;

                %e = abs(eigs(leftmat,rightmat,1,'SM','tolerance',1e-30));
                %e = abs(eigs(B*([M1+D1 sparse(n1,n2);sparse(n2,n1) M2+D2]\B'),Mlambda,1,'SM','tolerance',1e-30));
%                 betahs_h1h1 = [betahs_h1h1;sqrt(e)]
%                 figure(2)
%                 title('H1-H1')
%                 semilogy(i,e,'.','color',colors(count,:),'Markersize',10)
%                 if(e < 1e-16)
%                     break
%                 end
%                 hold on;
            end
        %catch e
        %    disp(e.message);
        %    break
        %end
    end
    
    betah_l2l2{end+1} = betahs_l2l2;
    betah_h1l2{end+1} = betahs_h1l2;
    betah_h1h1{end+1} = betahs_h1h1;
    nbfsh{end+1} = nbfs;
    condh{end+1} = conds;
end

if (exist('data','dir') == 0)
    disp('error')
    mkdir data;
end
%save('data/betah_h1h1.mat','betah_h1h1');
%save('data/betah_l2l2.mat','betah_l2l2');
%save('data/betah_h1l2.mat','betah_h1l2');
%save('data/nbfsh.mat','nbfsh');
%save('data/condh.mat','condh');
