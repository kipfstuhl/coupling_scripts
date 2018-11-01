clear all
close all
clc

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex')

alpha = 1;

% define exact solution and boundary conditions
u1ex = @(x) x(2,:).^2;
u2ex = @(x) cos(x(1,:));
pex = @(x) x(1,:).^2;

exsol_u = @(x) [u1ex(x);u2ex(x)];
exsol_p = @(x) pex(x);

u1exdx = @(x) 0;
u1exdy = @(x) 2*x(2,:).^1;
u1exdxdx = @(x) 0;
u1exdydy = @(x) 2;

u2exdx = @(x) -sin(x(1,:));
u2exdy = @(x) 0;
u2exdxdx = @(x) -cos(x(1,:));
u2exdydy = @(x) 0;

pexdx = @(x) 2*x(1,:);
pexdy = @(x) 0;

exsol_grad = @(x) [u1exdx(x) u1exdy(x);u2exdx(x) u2exdy(x)];

mu = 1;
fun = @(x) [-mu*(u1exdxdx(x)+u1exdydy(x)) + u1ex(x).*u1exdx(x) + u2ex(x).*u1exdy(x) + pexdx(x);
    -mu*(u2exdxdx(x)+u2exdydy(x)) + u1ex(x).*u2exdx(x) + u2ex(x).*u2exdy(x) + pexdy(x)];

dir = @(x) [exsol_u(x) exsol_u(x) exsol_u(x) exsol_u(x)  exsol_u(x)];
neu = @(x) zeros(2,6);

totalerr = [];
totalerr_lag = [];

for freqq = 0:1:5
    disp(['Freq = ',num2str(freqq)]);
    err = [];
    err_lag = [];
    erru = [];
    errp = [];
    hs = [];
    for ref = [1 3 5 7 9]
        ref_inflow = ref;
        ref_out1 = ref;
        ref_out2 = ref;
        
        mesh_in = read_mesh(['../meshes/refinement',num2str(ref_inflow),'/inflow_distorted.msh']);
        mesh_out1 = read_mesh(['../meshes/refinement',num2str(ref_out1),'/outflow1_distorted.msh']);
        mesh_out2 = read_mesh(['../meshes/refinement',num2str(ref_out2),'/outflow2_distorted.msh']);
        h = max(mesh_in.h,max(mesh_out1.h,mesh_out2.h));
        hs = [hs;h];
        fespace_us = cell(3,1);
        fespace_ps = cell(3,1);
        
        % fespaces inlet
        fespace_us{1} = create_fespace(mesh_in,'P2',[1 0 0 1 0]);
        fespace_ps{1} = create_fespace(mesh_in,'P1',[]);
        
        % fespaces outlet1
        fespace_us{2} = create_fespace(mesh_out1,'P2',[1 1 1 0 0]);
        fespace_ps{2} = create_fespace(mesh_out1,'P1',[]);
        
        % fespaces outlet2
        fespace_us{3} = create_fespace(mesh_out2,'P2',[0 1 1 1 0]);
        fespace_ps{3} = create_fespace(mesh_out2,'P1',[]);
        
        As = cell(3,1);
        bs = cell(3,1);
        
        [As{1},bs{1}] = assembler_steady_navier_stokes(fespace_us{1},fespace_ps{1},fun,mu,dir,neu);
        [As{2},bs{2}] = assembler_steady_navier_stokes(fespace_us{2},fespace_ps{2},fun,mu,dir,neu);
        [As{3},bs{3}] = assembler_steady_navier_stokes(fespace_us{3},fespace_ps{3},fun,mu,dir,neu);
        
        % find interface 1
        b1 = mesh_in.boundaries{2};
%         x1 = mesh_in.vertices(b1(1,1),1:2);
%         x2 = mesh_in.vertices(b1(end,end),1:2);
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

        nfreqs = [1 1 1] * freqq;

        n_nodes_us = cell(3,1);
        n_nodes_ps = cell(3,1);
        n_nodes_tot = cell(3,1);
        indices = cell(3,1);
        totalnodes = 0;
        cur = 0;
        for i = 1:3
            n_nodes_us{i} = size(fespace_us{i}.nodes,1);
            n_nodes_ps{i} = size(fespace_ps{i}.nodes,1);
            n_nodes_tot{i} = 2*n_nodes_us{i} + n_nodes_ps{i};
            indices{i} = cur+1:cur+2*n_nodes_us{i};
            cur = cur + n_nodes_tot{i};
            totalnodes = totalnodes + n_nodes_tot{i};
        end

        connectivity = [2 -5 0;
            0 4 -1;
            -3 0 5];
        
        n_boundaries = [5 5 5];
        
        xc = x1;
        
        B = [];
        
        prods = zeros(3);
        % add "constants"
        for i = 1:3 % interface 1
            bx = zeros(totalnodes,1);
            by = zeros(totalnodes,1);
            for j = 1:3 % interface 2
                prod = abs(ts(i,:)*ns(j,:)');
                prods(i,j) = prod;
                disp(['i = ', num2str(i),' j = ',num2str(j),' prod = ',num2str(prod)])
                for k = 1:3 % domains
                    flags = zeros(1,n_boundaries(k));
                    si = sign(connectivity(j,k));
                    if (si ~= 0)
                        flags(abs(connectivity(j,k))) = si;
                        b1 = zeros(n_nodes_us{k},1);
                        b2 = zeros(n_nodes_us{k},1);
                        b1 = apply_neumann_bc(b1,fespace_us{k},@(x) prod*flags,8);
                        bx(indices{k}) = bx(indices{k}) + [b1;b2];
                        by(indices{k}) = by(indices{k}) + [b2;b1];
                    end
                end
            end
            B = [B bx by];
        end
        
        for j = 1:3 % interface
            for i = 1:nfreqs(j)
                bx = zeros(totalnodes,2);
                by = zeros(totalnodes,2);
                
                f = @(x,xp,l) sin(i*sqrt((x(1,:)-xp(1)).^2+(x(2,:)-xp(2)).^2)/l*pi/2);
                for k = 1:3 % domains
                    flags = zeros(1,n_boundaries(k));
                    si = sign(connectivity(j,k));
                    if (si ~= 0)
                        flags(abs(connectivity(j,k))) = si;
                        b1 = zeros(n_nodes_us{k},1);
                        b2 = zeros(n_nodes_us{k},1);
                        b1 = apply_neumann_bc(b1,fespace_us{k},@(x) flags * f(x,xc,ls(j)),7);
                        bx(indices{k},1) = bx(indices{k},1) + [b1;b2];
                        by(indices{k},1) = by(indices{k},1) + [b2;b1];
                    end
                end
                
                f = @(x,xp,l) cos(i*sqrt((x(1,:)-xp(1)).^2+(x(2,:)-xp(2)).^2)/l*pi/2);
                for k = 1:3 % domains
                    flags = zeros(1,n_boundaries(k));
                    si = sign(connectivity(j,k));
                    if (si ~= 0)
                        flags(abs(connectivity(j,k))) = si;
                        b1 = zeros(n_nodes_us{k},1);
                        b2 = zeros(n_nodes_us{k},1);
                        b1 = apply_neumann_bc(b1,fespace_us{k},@(x) flags * f(x,xc,ls(j)),7);
                        bx(indices{k},2) = bx(indices{k},2) + [b1;b2];
                        by(indices{k},2) = by(indices{k},2) + [b2;b1];
                    end
                end
                B = [B bx by];
            end
        end
        B = sparse(B);
        
        Bt = B';
        dir_indices = [];
        cur = 0;
        for i = 1:3
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
        for i = 1:3
            A = @(u) blkdiag(A(u),As{i}(u(indices{i})));
            jac_block11 = @(u) blkdiag(jac_block11(u),build_jac_navier_stokes(As{i},u(indices{i}),fespace_us{i}));
        end
        
        x0 = zeros(length(bs{1}) + length(bs{2}) + length(bs{3}),1);
        x0 = [x0;zeros(nlamb,1)];
        
        jac = @(u) [jac_block11(u) B; Bt sparse(nlamb,nlamb)];
        
        H =@(u) [A(u) B;Bt sparse(nlamb,nlamb)];
        
        rhs = [bs{1};bs{2};bs{3};zeros(nlamb,1)];
        
        % solve system with newton's method
        ff = @(u) H(u)*u-rhs;
        tol = 1e-8;
        maxit = 20;
        
        [sol,er,it] = solve_with_newtons_method(ff,rhs,jac,tol,maxit);
        sols = cell(3,1);
        cur = 0;
        for i = 1:3
            sols{i}.u1 = sol(cur+1:cur+n_nodes_us{i});
            sols{i}.u2 = sol(cur+n_nodes_us{i}+1:cur+n_nodes_us{i}*2);
            sols{i}.p = sol(cur+2*n_nodes_us{i}+1:cur+n_nodes_tot{i});
            sols{i}.fespace_u = fespace_us{i};
            sols{i}.fespace_p = fespace_ps{i};
            cur = cur + n_nodes_tot{i};
        end
        
        lagrange_mult = sol(end-nlamb+1:end);
        npoints = 10000;
        prods_n = prods*lagrange_mult([1 3 5]);
        prods_t = prods*lagrange_mult([2 4 6]);
        errl = 0;
        cur = 6;
        for i = 1:3 % interface
            indices1_n = 1:4:nfreqs(i)*4;
            indices2_n = 2:4:nfreqs(i)*4;
            
            indices1_t = 3:4:nfreqs(i)*4;
            indices2_t = 4:4:nfreqs(i)*4;
            
            indices_n = reshape([indices1_n;indices2_n],nfreqs(i)*2,1)';
            indices_t = reshape([indices1_t;indices2_t],nfreqs(i)*2,1)';
            mylang_n = lagrange_mult(indices_n+cur);
            mylang_t = lagrange_mult(indices_t+cur);
            
            x = linspace(0,ls(i),npoints);
            st_exn_f = @(x) [1 0]*(mu*exsol_grad(x) * ns(i,:)' - exsol_p(x)*ns(i,:)');
            st_ext_f = @(x) [0 1]*(mu*exsol_grad(x) * ns(i,:)' - exsol_p(x)*ns(i,:)');
            
            st_exn = zeros(npoints,1);
            st_ext = zeros(npoints,1);
            for j = 1:npoints
                x_n = xc' + ts(i,:)'*x(j);
                st_exn(j) = st_exn_f(x_n);
                st_ext(j) = st_ext_f(x_n);
            end
            
            % reconstruct constant
            st_n = ones(npoints,1) * prods_n(i);
            st_t = ones(npoints,1) * prods_t(i);
            for j = 1:length(mylang_n)/2
                st_n = st_n + sin(j*x'/ls(i)*pi/2) * mylang_n((2*j-1));
                st_t = st_t + sin(j*x'/ls(i)*pi/2) * mylang_t((2*j-1));
                st_n = st_n + cos(j*x'/ls(i)*pi/2) * mylang_n(2*j);
                st_t = st_t + cos(j*x'/ls(i)*pi/2) * mylang_t(2*j);
            end
            
            dif_n = st_n - st_exn;
            dif_t = st_t - st_ext;
            
            errl = errl + sqrt(trapz(x,dif_n.^2)+trapz(x,dif_t.^2));
            if (nfreqs(i) ~= 0)
                cur = cur+indices2_t(end);
            end
        end
        
        err_lag = [err_lag errl]
        
        errsu = zeros(3,1);
        errsp = zeros(3,1);
        
        for i = 1:3
            errsu(i) = compute_H1_error_velocity(fespace_us{i},sols{i},exsol_u,exsol_grad);
            errsp(i) = compute_L2_error(fespace_ps{i},sols{i}.p,exsol_p);
        end
        
        errs = errsu + errsp;
        err = [err sqrt(errs'*errs)]
        erru = [erru sqrt(errsu'*errsu)];
        errp = [errp sqrt(errsp'*errsp)];
        
        vertices = [mesh_in.vertices;mesh_out1.vertices;mesh_out2.vertices];
        
        xm = min(vertices(:,1));
        xM = max(vertices(:,1));
        ym = min(vertices(:,2));
        yM = max(vertices(:,2));
        
        axis([xm xM ym yM])
        
        % plot_fe_function(sols{1}.u1,fespace_us{1})
        % plot_fe_function(sols{2}.u1,fespace_us{2})
        % plot_fe_function(sols{3}.u1,fespace_us{3})
        set(gca,'color','none')
        set(gca,'Visible','off')
        hold off
    end
    totalerr = [totalerr;err]
    totalerr_lag = [totalerr_lag;err_lag]
end

save('erroronlysincosdoublesupportlobatto.mat','totalerr')
save('erroronlysincosdoublesupportlobatto_lag.mat','totalerr_lag');

