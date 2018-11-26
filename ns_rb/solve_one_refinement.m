clear all
close all
clc

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex')

errsu = zeros(3,1);
errsp = zeros(3,1);

load('exsol.mat');
exsol = sol;

mu = 1;

exsol_u = @(x) [evaluate_fe_function(exsol.u1,exsol.fespace_u,x); ...
    evaluate_fe_function(exsol.u2,exsol.fespace_u,x)];
exsol_p = @(x) evaluate_fe_function(exsol.p,exsol.fespace_p,x);

exsol_grad = @(x) [evaluate_fe_function_gradient(exsol.u1,exsol.fespace_u,x)'; ...
    evaluate_fe_function_gradient(exsol.u2,exsol.fespace_u,x)'];

dir = @(x) [exsol_u(x) exsol_u(x) exsol_u(x) exsol_u(x)  exsol_u(x)];
neu = @(x) zeros(2,6);

totalerr = [];
totalerr_lag = [];

for freqq = 0
    disp(['Freq = ',num2str(freqq)]);
    err = [];
    err_lag = [];
    erru = [];
    errp = [];
    hs = [];
    for ref = 5
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
        fespace_us{1} = create_fespace(mesh_in,'P2',[1 0 0 1 1]);
        fespace_ps{1} = create_fespace(mesh_in,'P1',[]);
        
        % fespaces outlet1
        fespace_us{2} = create_fespace(mesh_out1,'P2',[1 0 1 0 0]);
        fespace_ps{2} = create_fespace(mesh_out1,'P1',[]);
        
        % fespaces outlet2
        fespace_us{3} = create_fespace(mesh_out2,'P2',[0 1 0 1 0]);
        fespace_ps{3} = create_fespace(mesh_out2,'P1',[]);
        
        r_in = 0.5;
        U = 1;
        dir0 = @(x) zeros(2,6);
        dir_in = @(x) [0 0; 0 0; 0 0;0 0; (r_in^2 - x(2)^2)/r_in^2 * U 0]';
        neu0 = dir0;
        
        
        As = cell(3,1);
        bs = cell(3,1);
        
        [As{1},bs{1}] = assembler_steady_navier_stokes(fespace_us{1},fespace_ps{1},[0;0],1,dir_in,neu0);
        [As{2},bs{2}] = assembler_steady_navier_stokes(fespace_us{2},fespace_ps{2},[0;0],1,dir0,neu0);
        [As{3},bs{3}] = assembler_steady_navier_stokes(fespace_us{3},fespace_ps{3},[0;0],1,dir0,neu0);
        
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
        %                 for i = 1:3 % interface 1
        %                     bx = zeros(totalnodes,1);
        %                     by = zeros(totalnodes,1);
        %                     for j = 1:3 % interface 2
        %                         prod = abs(ts(i,:)*ns(j,:)');
        %                         prods(i,j) = prod;
        %                         disp(['i = ', num2str(i),' j = ',num2str(j),' prod = ',num2str(prod)])
        %                         for k = 1:3 % domains
        %                             flags = zeros(1,n_boundaries(k));
        %                             si = sign(connectivity(j,k));
        %                             if (si ~= 0)
        %                                 flags(abs(connectivity(j,k))) = si;
        %                                 b1 = zeros(n_nodes_us{k},1);
        %                                 b2 = zeros(n_nodes_us{k},1);
        %                                 b1 = apply_neumann_bc(b1,fespace_us{k},@(x) prod*flags,8);
        %                                 bx(indices{k}) = bx(indices{k}) + [b1;b2];
        %                                 by(indices{k}) = by(indices{k}) + [b2;b1];
        %                             end
        %                         end
        %                     end
        %                     B = [B bx by];
        %                 end
        
        for j = 1:3 % interface
            for i = 0:nfreqs(j)
                if (i == 0)
                    bx = zeros(totalnodes,1);
                    by = zeros(totalnodes,1);
                    f = @(x,xp,l) 1;
                    for k = 1:3 % domains
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
                    for k = 1:3 % domains
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
                    for k = 1:3 % domains
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
        
        hold on
        unorm1 = sqrt(sols{1}.u1.^2 + sols{1}.u2.^2);
        unorm2 = sqrt(sols{2}.u1.^2 + sols{2}.u2.^2);
        unorm3 = sqrt(sols{3}.u1.^2 + sols{3}.u2.^2);
        
        minnorm = min([unorm1;unorm2;unorm3]);
        maxnorm = max([unorm1;unorm2;unorm3]);
        
        plot_fe_fluid_function(sols{1},'U')
        plot_fe_fluid_function(sols{2},'U')
        plot_fe_fluid_function(sols{3},'U')
        caxis([minnorm maxnorm]);
        
        vertices = [mesh_in.vertices;mesh_out1.vertices;mesh_out2.vertices];
        
        xm = min(vertices(:,1));
        xM = max(vertices(:,1));
        ym = min(vertices(:,2));
        yM = max(vertices(:,2));
        
        axis([xm xM ym yM])
        
        axis equal
        set(gca,'color','none')
        set(gca,'Visible','off')
        hold off
        colorbar('Location','south')
        
        figure
        
        load(['solutions/u1s_ref',num2str(ref),'.mat']);
        load(['solutions/u2s_ref',num2str(ref),'.mat']);
        load(['solutions/ps_ref',num2str(ref),'.mat']);
        
        bcmatrix = [1 0 0 1 1; 1 1 1 0 0; 0 1 1 1 0];
        minlog = Inf;
        maxlog = -Inf;
        for i = 1:3
            fespace_us{i}.bc = bcmatrix(i,:);
            
            indu = find_dirichlet_indices(fespace_us{i});
            
            normu = sqrt(sols{i}.u1.^2 + sols{i}.u2.^2);
            
            dif1 = abs(sols{i}.u1 - u1s{i});
            dif1(indu) = 1e-6;
            
            dif2 = abs(sols{i}.u2 - u2s{i});
            dif2(indu) = 1e-6;

            N = log(sqrt(dif1.^2 + dif2.^2));
           
            hold on
            [~,h] = tricontf(fespace_us{i}.mesh.vertices(:,1),fespace_us{i}.mesh.vertices(:,2), ...
                   fespace_us{i}.mesh.elements(:,1:3),N(1:size(fespace_us{i}.mesh.vertices(:,1),1)),20);
            set(h,'edgecolor','none');
            minlog = min(min(N),minlog);
            maxlog = max(max(N),maxlog);
        end
        c = linspace(minlog,maxlog,6);
        expc = exp(c);
        for i = 1:length(c)
            check = 1;
            count = 0;
            while check
                count = count + 1;
                if (floor(expc(i) * 10^count) > 0)
                    check = false;
                end
            end
            expc(i) = 10^(-count) * floor(expc(i) * 10^count);
        end
        
        colorbar('YTick',c,'YTickLabel',expc,'Location','south');
        axis([xm xM ym yM])
        
        axis equal
        set(gca,'color','none')
        set(gca,'Visible','off')
        hold off
    end
    totalerr = [totalerr;err]
    totalerr_lag = [totalerr_lag;err_lag]
end

