clc
clear all
close all

run load_exact_solution_poisson.m

dir_functions = @(x) [0;0;0;0];
neu_functions = @(x) [0;0;0;0];

N = [20 28 40 56 80 114 160];

nrefs = length(N);
nfreq = 15;
condnum = zeros(nfreq,nrefs+1);

% compute matrix for the orthonormalization of the fourier basis f
R = compute_orthonormalization_matrix(nfreq,10000);

beta = {};
count = 0;
for n_elements = N
    count = count + 1;
    n1x = n_elements/2;
    n2x = n_elements/2;
    n1y = n_elements;
    n2y = n_elements+1;
    
    xp1 = 0;
    yp1 = 0;
    L1 = 0.5;
    H1 = 1;
    
    % create mesh 1
    mesh1 = create_mesh(xp1,yp1,L1,H1,n1x,n1y);
    
    xp2 = L1;
    yp2 = 0;
    L2 = 0.5;
    H2 = 1;
    
    % create mesh 2
    mesh2 = create_mesh(xp2,yp2,L2,H2,n2x,n2y);
    
    % create fespace 1
    bc1 = [1 0 1 1];
    fespace1 = create_fespace(mesh1,'P2',bc1);
    
    % create interface mass for domain 1
    [M1,RS1] = build_interface_matrix_and_restriction(fespace1,2);
    
    % create fespace 2
    bc2 = [1 1 1 0];
    fespace2 = create_fespace(mesh2,'P2',bc2);
    
    % create interface mass for domain 2
    [M2,RS2] = build_interface_matrix_and_restriction(fespace2,4);
    
    n1 = size(M1,1);
    n2 = size(M2,1);
        
    B1 = [];
    B2 = [];
    
    betas = [];
    for i = 0:nfreq
        disp(['N elements = ',num2str(n_elements),', freq = ',num2str(i)]);
        B1 = add_row_to_coupling_matrix_poisson(B1,fespace1,2,i,4);
        B2 = add_row_to_coupling_matrix_poisson(B2,fespace2,4,i,4);
        
        nlm = size(B1,1);
        
        % restrict R to the basis functions that we are actually using
        Rpar = R(1:nlm,1:nlm);
        
        B1GS = Rpar'\B1*RS1';
        B2GS = Rpar'\B2*RS2';
        
        % assemble matrices of the eigenvalue problem
        M = [M1 sparse(n1,n2);
               sparse(n2,n1) M2];
        
        RB = [-(Rpar'\B1)*RS1' (Rpar'\B2)*RS2'];   
        mat = full(RB*(M\RB'));

        % compute minimum eigenvalue and its sqrt
        e = sqrt(min(eig(mat)));
        ngamma = 2*i+1;

        betas = [betas;[ngamma e]];
    end
    beta{end+1} = betas;
end

% save data
save('data_figure5/beta.mat','beta');
