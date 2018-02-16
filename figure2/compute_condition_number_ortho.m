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

R = compute_orthonormalization_matrix(nfreq,10000);

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
    
    [A1,rhs1] = assembler_poisson(fespace1,fun,1,dir_functions,neu_functions);
    
    bc2 = [1 1 1 0];
    fespace2 = create_fespace(mesh2,'P2',bc2);
    
    [A2,rhs2] = assembler_poisson(fespace2,fun,1,dir_functions,neu_functions);
    
    n1 = size(A1,1);
    n2 = size(A2,1);
    
    indices1 = 1:n1;
    indices2 = n1+1:n1+n2;
        
    B1 = [];
    B2 = [];
    
    for i = 0:nfreq
        disp(['N elements = ',num2str(n_elements),', freq = ',num2str(i)]);
        B1 = add_row_to_coupling_matrix_poisson(B1,fespace1,2,i);
        B2 = add_row_to_coupling_matrix_poisson(B2,fespace2,4,i);
        
        nlm = size(B1,1);
        
        % restrict R to the basis functions that we are actually using
        Rpar = R(1:nlm,1:nlm);
        
        B1GS = Rpar'\B1;
        B2GS = Rpar'\B2;
        
        mat = [A1 sparse(n1,n2) -B1GS';
               sparse(n2,n1) A2  B2GS';
               -B1GS B2GS sparse(nlm,nlm)];
           
        k = condest(mat);
        condnum(i+1,count) = k;
    end
end

save('data_figure2/condnum_ortho.mat','condnum');
