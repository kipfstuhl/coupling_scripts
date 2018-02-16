function beta = compute_beta(poly1,poly2,norm1,norm2)

run load_exact_solution_and_f.m

dir_functions = @(x) [0;0;0;0];
neu_functions = @(x) [0;0;0;0];

N = [20 28 40 56 80 114 160];
% N = [20 28 40 56 80];

colors = [1 0 0;
    0 1 0;
    0 0 1;
    0.5 0.5 0;
    0 0.5 0.5];

beta = {};

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
    fespace1 = create_fespace(mesh1,poly1,bc1);
    
    [A1,rhs1] = assembler_poisson(fespace1,fun,1,dir_functions,neu_functions);
    
    % create also mass matrix for first problem (needed to compute norm
    % matrix)
    %M1 = apply_dirichlet_bc_matrix(assemble_mass(fespace1),fespace1,1);
    D1 = assemble_stiffness(1,fespace1);
    M1 = assemble_mass(fespace1);
    
    bc2 = [1 1 1 0];
    fespace2 = create_fespace(mesh2,poly2,bc2);
    
    [A2,rhs2] = assembler_poisson(fespace2,fun,1,dir_functions,neu_functions);
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
    
    betas = [];
    bfs = {};
    bfsd = {};
    
    
    for i = 1:n_iterations
        disp(['n_elements = ', num2str(n_elements) ...
            ', iteration n = ', num2str(i)]);
        
        b1 = zeros(n1,1);
        b2 = zeros(n2,1);
        
        freq = (i - 1)*2;
        coef = sqrt(2);
        
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
            b1 = apply_neumann_bc(b1,fespace1,@(x) [0;coef*sin(x(2) * pi * freq);0;0],quadrature_order);
            B1 = [B1;b1'];
            
            B1_t = [B1_t apply_dirichlet_bc_rhs(b1,fespace1,@(x) [0;0;0;0])];
            
            
            b1 = b1*0;
            
            bfs{end+1} = @(x) coef*sin(x * pi * freq);
            bfsd{end+1} = @(x) coef*pi * freq*cos(x * pi * freq);
            
            b1 = apply_neumann_bc(b1,fespace1,@(x) [0;coef*cos(x(2) * pi * freq);0;0],quadrature_order);
            B1 = [B1;b1'];
            
            B1_t = [B1_t apply_dirichlet_bc_rhs(b1,fespace1,@(x) [0;0;0;0])];
            
            
            bfs{end+1} = @(x) coef*cos(x * pi * freq);
            bfsd{end+1} = @(x) -coef*pi * freq*sin(x * pi * freq);
            
            b2 = apply_neumann_bc(b2,fespace2,@(x) [0;0;0;coef*sin(x(2) * pi * freq)],quadrature_order);
            B2 = [B2;b2'];
            
            B2_t = [B2_t apply_dirichlet_bc_rhs(b2,fespace2,@(x) [0;0;0;0])];
            
            b2 = b2*0;
            
            b2 = apply_neumann_bc(b2,fespace2,@(x) [0;0;0;coef*cos(x(2) * pi * freq)],quadrature_order);
            B2 = [B2;b2'];
            
            B2_t = [B2_t apply_dirichlet_bc_rhs(b2,fespace2,@(x) [0;0;0;0])];
        end
        
        n = n1 + n2;
        nm = size(B1,1);
        
        %try
            lambdas = size(bfs,2);
            Mlambda = speye(lambdas);
            
            B = [-B1 B2];
                        
            if (strcmp(norm2,'H1'))
                error('Norm H1 for the Lagrange multiplier is not implemented');
            end
            mat1 = [M1+D1 sparse(n1,n2) -B1'; sparse(n2,n1) M2+D2 B2'; -B1 B2 sparse(lambdas,lambdas)];
            mat2 = [sparse(n1,n1) sparse(n1,n2) sparse(n1,lambdas);
                    sparse(n2,n1) sparse(n2,n2) sparse(n2,lambdas);
                    sparse(lambdas,n1) sparse(lambdas,n2) Mlambda];
            %mat = full(B*([M1+D1 sparse(n1,n2); sparse(n2,n1) M2+D2]\B'));
            e = -eigs(mat1,mat2,1,'SM');
            % e = min(eig(mat));
            if (e < 1e-16)
                break
            end
            e = sqrt(e);
            ngamma = 2*i-1;
            %             figure(1)
            %             semilogy(ngamma,sqrt(e),'.','color',colors(count,:),'Markersize',10)
            %             hold on
            %             xlim([0 ngamma])
            
            betas = [betas;[ngamma sqrt(e)]];
%             
%         catch e
%             disp(e.message);
%             break
%         end
    end
    beta{end+1} = betas;
end
