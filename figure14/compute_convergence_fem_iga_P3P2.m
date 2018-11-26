clear all
close all
clc

load('../fine_solution_bifurcation/exsol.mat');
u1ex = @(x,y) evaluate_fe_function(exsol.u1,exsol.fespace_u,[x(:)';y(:)']);
u2ex = @(x,y) evaluate_fe_function(exsol.u2,exsol.fespace_u,[x(:)';y(:)']);

% load interpolated solutions for
load('data_figure14/u1s_ref2.mat');
load('data_figure14/u2s_ref2.mat');
load('data_figure14//ps_ref2.mat');

velex = @(x, y) cat(1, ...
    reshape (u1ex (x,y), [1, size(x)]), ...
    reshape (u2ex (x,y), [1, size(x)]));

gradex = @(x,y) cat(1, ...
    reshape (evaluate_fe_function_gradient(exsol.u1,exsol.fespace_u,[x(:)';y(:)']),[1,2,size(x)]), ...
    reshape (evaluate_fe_function_gradient(exsol.u2,exsol.fespace_u,[x(:)';y(:)']),[1,2,size(x)]));

pex = @(x,y) evaluate_fe_function(exsol.p,exsol.fespace_p,[x(:)';y(:)']);

% degree of polynomials for pressure to be used in outflow 
deg = 2;

% options inflow (IGA)
options_inflow.element_name = 'th';     % Element type for discretization
options_inflow.degree = [ 3  3];        % Degree of the splines (pressure space)
options_inflow.regularity = [ 2  2];    % Regularity of the splines (pressure space)
options_inflow.nsub = [ 6  6];          % Number of subdivisions
options_inflow.nquad = [ 5  5];         % Points for the Gaussian quadrature rule

% options outflow 1 (FEM)
options_outflow1.polydegree_u = 'P2';
options_outflow1.polydegree_p = 'P1';
options_outflow1.ref = 2;

count = 0;
for nel_iga = [5:2:15]
    count = count + 1;
    for freq = 0:5
        
        % options outflow 2 (IGA)
        options_outflow2.element_name = 'th';
        options_outflow2.degree = [ deg  deg];
        options_outflow2.regularity = [ deg-1  deg-1];
        options_outflow2.nsub = [ nel_iga  nel_iga];
        options_outflow2.nquad = [ 5  5];
        
        [sol_fem,is] = solve_system_bifurcation_FEM_IGA(options_inflow,options_outflow1,options_outflow2,freq);
       
        % compute error on the inflow 
        errsp = sp_l2_error (is{1}.space_p, is{1}.msh, is{1}.p, pex);
        errsu = sp_h1_error (is{1}.space_v, is{1}.msh, is{1}.vel, velex, gradex);
        
        err1(count,freq+1) = errsu + errsp
        
        % compute finite element error against interpolated fine function      
        dif1 = sol_fem.u1 - u1s{2};
        e1 = compute_H1_error(sol_fem.fespace_u,dif1,@(x) 0,@(x) [0;0]);
        
        dif2 = sol_fem.u2 - u2s{2};
        e2 = compute_H1_error(sol_fem.fespace_u,dif2,@(x) 0,@(x) [0;0]);

        errsu = sqrt(e1^2 + e2^2); 
  
        difp = sol_fem.p - ps{2};
        errsp = compute_L2_error(sol_fem.fespace_p,difp,@(x) 0);
          
        err2(count,freq+1) = errsu + errsp
        
        % retrieving mesh size
        msh_prc = msh_precompute(is{2}.msh);
        h_iga(count) = max(msh_prc.element_size);
    end
end

save(['data_figure14/errs_P3P2'],'err1','err2','h_iga')

