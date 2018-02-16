function cond = compute_conds(poly1,poly2)

run load_exact_solution_and_f.m

dir_functions = @(x) [0;0;0;0];
neu_functions = @(x) [0;0;0;0];

N = [20 28 40 56 80];

cond = {};

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
    
    [A1,rhs1] = assembler_poisson(fespace1,fun,mu,dir_functions,neu_functions);
    
    % create also mass matrix for first problem (needed to compute norm
    % matrix)
    %M1 = apply_dirichlet_bc_matrix(assemble_mass(fespace1),fespace1,1);
    D1 = assemble_stiffness(1,fespace1);
    M1 = assemble_mass(fespace1);
    
    bc2 = [1 1 1 0];
    fespace2 = create_fespace(mesh2,poly2,bc2);
    
    [A2,rhs2] = assembler_poisson(fespace2,fun,mu,dir_functions,neu_functions);
    %M2 = apply_dirichlet_bc_matrix(assemble_mass(fespace2),fespace2,1);
    D2 = assemble_stiffness(1,fespace2);
    M2 = assemble_mass(fespace2);
    
    n1 = size(A1,1);
    n2 = size(A2,1);
    
    indices1 = 1:n1;
    indices2 = n1+1:n1+n2;
    
    n_iterations = 16;
    
    [basis,xx] = create_ortho_fourier_basis(n_iterations,1000);
    nbasis = size(basis,2);
    
    evaluate_basis = @(x,vec) interp1(xx,vec,x);
    
    B1 = [];
    B2 = [];
    
    B1_t = [];
    B2_t = [];
    
    conds = [];
    bfs = {};
    bfsd = {};
    
    for i = 1:nbasis
        disp(['n_elements = ', num2str(n_elements) ...
            ', n basis functions = ', num2str(i)]);
        
        b1 = zeros(n1,1);
        b2 = zeros(n2,1);
        
        b1 = apply_neumann_bc(b1,fespace1,@(x) [0;evaluate_basis(x(2),basis(:,i));0;0]);

        B1 = [B1;b1'];

        B1_t = [B1_t apply_dirichlet_bc_rhs(b1,fespace1,@(x) [0;0;0;0])];

        b2 = apply_neumann_bc(b2,fespace2,@(x) [0;0;0;evaluate_basis(x(2),basis(:,i));]);

        B2 = [B2;b2'];
        B2_t = [B2_t apply_dirichlet_bc_rhs(b2,fespace2,@(x) [0;0;0;0])];
        
        n = n1 + n2;
        nm = size(B1,1);
        
        try
           
            
            mat = [A1 sparse(n1,n2) -B1_t;
                sparse(n2,n1) A2 B2_t;
                -B1 B2 sparse(nm,nm)];
            
            ngamma = i;
            
            if (mod(ngamma,2) == 1)
                conds = [conds;[ngamma condest(mat)]];
            end
            if (conds(end,2) > 1e20)
                break;
            end
            
        catch e
            disp(e.message);
            break
        end
    end
    cond{end+1} = conds;
end
