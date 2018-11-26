clear all
close all
clc

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex')

errsu = zeros(3,1);
errsp = zeros(3,1);

mu = 1;

dir = @(x) [exsol_u(x) exsol_u(x) exsol_u(x) exsol_u(x)  exsol_u(x)];
neu = @(x) zeros(2,6);

totalerr = [];
totalerr_lag = [];

ref_inflow = 12;
ref_out1 = 12;

countdeg = 0;

% degree of polynomials to be used in outflow 2
deg = 2;

count = 0;
for nel_iga = [5:2:15]
    count = count + 1;
    for freq = 0:5
        % ========================= FINITE ELEMENT PART =========================
        
        % load finite element matrices (we load the second outflow just to retrieve
        % geometric information about the interface 3)
        mesh_in = read_mesh(['../meshes/refinement',num2str(ref_inflow),'/inflow_distorted.msh']);
        mesh_out1 = read_mesh(['../meshes/refinement',num2str(ref_out1),'/outflow1_distorted.msh']);
        mesh_out2 = read_mesh(['../meshes/refinement',num2str(ref_out1),'/outflow2_distorted.msh']);
        h = max(mesh_in.h,mesh_out1.h);
        
        fespace_us = cell(2,1);
        fespace_ps = cell(2,1);
        
        % fespaces outlet1
        fespace_us{2} = create_fespace(mesh_out1,'P2',[1 0 1 0 0]);
        fespace_ps{2} = create_fespace(mesh_out1,'P1',[]);
        
        r_in = 0.5;
        U = 1;
        dir0 = @(x) zeros(2,6);
        dir_in = @(x) [0 0; 0 0; 0 0;0 0; (r_in^2 - x(2)^2)/r_in^2 * U 0]';
        neu0 = dir0;
        
        As = cell(3,1);
        bs = cell(3,1);
        
        [As{2},bs{2}] = assembler_steady_navier_stokes(fespace_us{2},fespace_ps{2},[0;0],1,dir0,neu0);
        
        % find interface 1
        b1 = mesh_in.boundaries{2};
        x1 = mesh_in.vertices(b1(end,end),1:2);
        x2 = mesh_in.vertices(b1(1,1),1:2);
        dif = x2 - x1;
        l1 = norm(dif);
        t1 = dif/l1;
        n1 = [t1(2) -t1(1)];
        
        % find interface 2
        b2 = mesh_out1.boundaries{4};
        x1 = mesh_out1.vertices(b2(1,1),1:2);
        x2 = mesh_out1.vertices(b2(end,end),1:2);
        dif = x2 - x1;
        l2 = norm(dif);
        t2= dif/l2;
        n2 = [t2(2) -t2(1)];
        
        % find interface 3
        b3 = mesh_out2.boundaries{5};
        x1 = mesh_out2.vertices(b3(1,1),1:2);
        x2 = mesh_out2.vertices(b3(end,end),1:2);
        dif = x2 - x1;
        l3 = norm(dif);
        t3= dif/l3;
        n3 = [t3(2) -t3(1)];
        
        ts = [t1;t2;t3];
        ns = [n1;n2;n3];
        ls = [l1 l2 l3];
        
        nfreqs = [1 1 1] * freq;
        
        n_nodes_us = cell(2,1);
        n_nodes_ps = cell(2,1);
        n_nodes_tot = cell(2,1);
        indices = cell(2,1);
        
        totalnodes = 0;
        cur = 0;
        for i = 2:2
            n_nodes_us{i} = size(fespace_us{i}.nodes,1);
            n_nodes_ps{i} = size(fespace_ps{i}.nodes,1);
            n_nodes_tot{i} = 2*n_nodes_us{i} + n_nodes_ps{i};
            indices{i} = cur+1:cur+2*n_nodes_us{i};
            cur = cur + n_nodes_tot{i};
            totalnodes = totalnodes + n_nodes_tot{i};
        end
        
        connectivity = [2 -5;
            0 -4;
            -3 0];
        
        n_boundaries = [5 5];
        
        xc = x1;
        
        B = [];
        
        prods = zeros(3);
        for j = 1:3 % interface
            for i = 0:nfreqs(j)
                if (i == 0)
                    bx = zeros(totalnodes,1);
                    by = zeros(totalnodes,1);
                    f = @(x,xp,l) 1;
                    for k = 2:2 % domains
                        flags = zeros(1,n_boundaries(k));
                        si = sign(connectivity(j,k));
                        if (si ~= 0)
                            flags(abs(connectivity(j,k))) = si;
                            b1 = zeros(n_nodes_us{k},1);
                            b2 = zeros(n_nodes_us{k},1);
                            b1 = apply_neumann_bc(b1,fespace_us{k},@(x) flags * f(x,xc,ls(j)),8);
                            bx(indices{k},1) = bx(indices{k},1) + [b1;b2];
                            by(indices{k},1) = by(indices{k},1) + [b2;b1];
                        end
                    end
                else
                    bx = zeros(totalnodes,2);
                    by = zeros(totalnodes,2);
                    f = @(x,xp,l) sin(i*sqrt((x(1,:)-xp(1)).^2+(x(2,:)-xp(2)).^2)/l*pi);
                    for k = 2:2 % domains
                        flags = zeros(1,n_boundaries(k));
                        si = sign(connectivity(j,k));
                        if (si ~= 0)
                            flags(abs(connectivity(j,k))) = si;
                            b1 = zeros(n_nodes_us{k},1);
                            b2 = zeros(n_nodes_us{k},1);
                            b1 = apply_neumann_bc(b1,fespace_us{k},@(x) flags * f(x,xc,ls(j)),8);
                            bx(indices{k},1) = bx(indices{k},1) + [b1;b2];
                            by(indices{k},1) = by(indices{k},1) + [b2;b1];
                        end
                    end
                    
                    f = @(x,xp,l) cos(i*sqrt((x(1,:)-xp(1)).^2+(x(2,:)-xp(2)).^2)/l*pi);
                    for k = 2:2 % domains
                        flags = zeros(1,n_boundaries(k));
                        si = sign(connectivity(j,k));
                        if (si ~= 0)
                            flags(abs(connectivity(j,k))) = si;
                            b1 = zeros(n_nodes_us{k},1);
                            b2 = zeros(n_nodes_us{k},1);
                            b1 = apply_neumann_bc(b1,fespace_us{k},@(x) flags * f(x,xc,ls(j)),8);
                            bx(indices{k},2) = bx(indices{k},2) + [b1;b2];
                            by(indices{k},2) = by(indices{k},2) + [b2;b1];
                        end
                    end
                end
                B = [B bx by];
            end
        end
        B = sparse(B);
        
        Bt = B';
        dir_indices = [];
        cur = 0;
        for i = 2:2
            loc_dir_indices = find_dirichlet_indices(fespace_us{i});
            dir_indices = [dir_indices;
                loc_dir_indices + cur;
                loc_dir_indices + cur + n_nodes_us{i}];
            cur = cur + n_nodes_tot{i};
        end
        B(dir_indices,:) = 0;
        
        nlamb = size(B,2);
        
        jac_block11 = @(u) [];
        A = @(u) [];
        for i = 2:2
            A = @(u) blkdiag(A(u),As{i}(u(indices{i})));
            jac_block11 = @(u) blkdiag(jac_block11(u),build_jac_navier_stokes(As{i},u(indices{i}),fespace_us{i}));
        end
        
        A_iga = cell(2,1);
        B_iga = cell(2,1);
        J_iga = cell(2,1);
        rhs_iga = cell(2,1);
        nurb = cell(2,1);
        int_dofs = cell(2,1);
        vel_iga = cell(2,1);
        Ndofs = cell(2,1);
        space_v = cell(2,1);
        space_p = cell(2,1);
        geometry = cell(2,1);
        msh = cell(2,1);
        dofs_iga = cell(2,1);
        % ========================= IGA PART INFLOW =========================
        
        method_data.element_name = 'th';     % Element type for discretization
        method_data.degree       = [ 3  3];  % Degree of the splines (pressure space)
        method_data.regularity   = [ 2  2];  % Regularity of the splines (pressure space)
        method_data.nsub         = [ 40  40];  % Number of subdivisions
        method_data.nquad        = [ 5  5];  % Points for the Gaussian quadrature rule
        
        ii = 1;
        [A_iga{ii},B_iga{ii},J_iga{ii},rhs_iga{ii},nurb{ii},int_dofs{ii},vel_iga{ii},Ndofs{ii},space_v{ii},space_p{ii},geometry{ii},msh{ii}] = ...
            assemble_system_iga_inflow(method_data,ls,nfreqs);
        keyboard
        % ========================= IGA PART OUTFLOW 2 =========================
        
        method_data.element_name = 'th';     % Element type for discretization
        method_data.degree       = [ deg  deg];  % Degree of the splines (pressure space)
        method_data.regularity   = [ deg-1  deg-1];  % Regularity of the splines (pressure space)
        method_data.nsub         = [ nel_iga  nel_iga];  % Number of subdivisions
        method_data.nquad        = [ 5  5];  % Points for the Gaussian quadrature rule
        
        % Penalization parameter for Nitsche's method
        factor = 10;
        method_data.Cpen = factor*(min(method_data.degree)+1);
        
        ii = 2;
        [A_iga{ii},B_iga{ii},J_iga{ii},rhs_iga{ii},nurb{ii},int_dofs{ii},vel_iga{ii},Ndofs{ii},space_v{ii},space_p{ii},geometry{ii},msh{ii}] = ...
            assemble_system_iga_outflow2(method_data,ls,nfreqs);
        
        aux = totalnodes;
        for i = 1:2
            B = [B;B_iga{i}];
            Bt = [Bt B_iga{i}'];
            dofs_iga{i} = aux+1:aux+Ndofs{i};
            A = @(u) blkdiag(A(u),A_iga{i}(u(dofs_iga{i})));
            jac_block11 = @(u) blkdiag(jac_block11(u),J_iga{i}(u(dofs_iga{i})));
            aux = aux+Ndofs{i};
        end
        
        H =@(u) [A(u) B;Bt sparse(nlamb,nlamb)];
        jac = @(u) [jac_block11(u) B; Bt sparse(nlamb,nlamb)];
        
        x0 = zeros(length(bs{2}) + Ndofs{1} + Ndofs{2} + nlamb,1);
        
        rhs = @(u) [bs{2};rhs_iga{1}(u(dofs_iga{1}));rhs_iga{2}(u(dofs_iga{2}));zeros(nlamb,1)];
        
        % solve system with newton's method
        ff = @(u) H(u)*u-rhs(u);
        tol = 1e-8;
        maxit = 20;
        
        [sol,er,it] = solve_with_newtons_method(ff,x0,jac,tol,maxit);
        
        sols_fem = cell(2,1);
        cur = 0;
        for i = 2:2
            sols_fem{i}.u1 = sol(cur+1:cur+n_nodes_us{i});
            sols_fem{i}.u2 = sol(cur+n_nodes_us{i}+1:cur+n_nodes_us{i}*2);
            sols_fem{i}.p = sol(cur+2*n_nodes_us{i}+1:cur+n_nodes_tot{i});
            sols_fem{i}.fespace_u = fespace_us{i};
            sols_fem{i}.fespace_p = fespace_ps{i};
            cur = cur + n_nodes_tot{i};
        end
        
        sol_iga = sol(dofs_iga{1});
        v_dofs = sol_iga(1:length(int_dofs{1}));
        vel_iga{1}(int_dofs{1}) = v_dofs;
        p_iga{1} = sol_iga(length(int_dofs{1})+1:end);
        
        sol_iga = sol(dofs_iga{2});
        v_dofs = sol_iga(1:length(int_dofs{2}));
        vel_iga{2}(int_dofs{2}) = v_dofs;
        p_iga{2} = sol_iga(length(int_dofs{2})+1:end);
        
        dofsu1 = 2 * numel(sols_fem{2}.u1);
        dofsu2 = numel(vel_iga{1});
        dofsu3(count) = numel(vel_iga{2});
        
        dofsp1 = 2 * numel(sols_fem{2}.p);
        dofsp2 = numel(p_iga{1});
        dofsp3(count) = numel(p_iga{2});
        
        % ========================= COMPUTING ERRORS =========================
        
        disp('Outflow 1');
        % outflow 1
        
        load(['../ns_rb/solutions/u1s_ref',num2str(ref_out1),'.mat']);
        load(['../ns_rb/solutions/u2s_ref',num2str(ref_out1),'.mat']);
        load(['../ns_rb/solutions/ps_ref',num2str(ref_out1),'.mat']);
        
        bcmatrix = [1 0 0 1 1; 1 1 1 0 0; 0 1 1 1 0];
        for i = 2
            dif1 = sols_fem{i}.u1 - u1s{i};
            e1 = compute_H1_error(fespace_us{i},dif1,@(x) 0,@(x) [0;0]);
            dif2 = sols_fem{i}.u2 - u2s{i};
            e2 = compute_H1_error(fespace_us{i},dif2,@(x) 0,@(x) [0;0]);
            errsu = sqrt(e1^2 + e2^2);       
            difp = sols_fem{i}.p - ps{i};
            e = compute_L2_error(fespace_ps{i},difp,@(x) 0);
            errsp = e;
        end
        
        err1(count,freq+1) = errsu + errsp
        
        load('../ns_rb/solutions/exsol.mat');
        exsol = sol;
        
        u1ex = @(x,y) evaluate_fe_function(exsol.u1,exsol.fespace_u,[x(:)';y(:)']);
        u2ex = @(x,y) evaluate_fe_function(exsol.u2,exsol.fespace_u,[x(:)';y(:)']);
        
        velex = @(x, y) cat(1, ...
            reshape (u1ex (x,y), [1, size(x)]), ...
            reshape (u2ex (x,y), [1, size(x)]));
        
        gradex = @(x,y) cat(1, ...
            reshape (evaluate_fe_function_gradient(exsol.u1,exsol.fespace_u,[x(:)';y(:)']),[1,2,size(x)]), ...
            reshape (evaluate_fe_function_gradient(exsol.u2,exsol.fespace_u,[x(:)';y(:)']),[1,2,size(x)]));
        
        pex = @(x,y) evaluate_fe_function(exsol.p,exsol.fespace_p,[x(:)';y(:)']);
        
        disp('Inflow');
        % inflow 
        errsp = sp_l2_error (space_p{1}, msh{1}, p_iga{1}, pex);
        errsu = sp_h1_error (space_v{1}, msh{1}, vel_iga{1}, velex, gradex);
        
        err2(count,freq+1) = errsu + errsp
        
        disp('Outflow 2');
        % outflow2 
        errsp = sp_l2_error (space_p{2}, msh{2}, p_iga{2}, pex);
        errsu = sp_h1_error (space_v{2}, msh{2}, vel_iga{2}, velex, gradex);
        
        err3(count,freq+1) = errsu + errsp 
        
        msh_prc = msh_precompute(msh{2});
        h_iga(count) = max(msh_prc.element_size)
    end
end

save(['errs_',method_data.element_name],'_P3P2','dofsu1','dofsu2','dofsu3', ...
    'dofsp1','dofsp2','dofsp3','err1','err2','err3')

