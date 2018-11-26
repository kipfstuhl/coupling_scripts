clear all
close all
clc

load('../fine_solution_bifurcation/exsol.mat');
u1ex = @(x,y) evaluate_fe_function(exsol.u1,exsol.fespace_u,[x(:)';y(:)']);
u2ex = @(x,y) evaluate_fe_function(exsol.u2,exsol.fespace_u,[x(:)';y(:)']);

velex = @(x, y) cat(1, ...
    reshape (u1ex (x,y), [1, size(x)]), ...
    reshape (u2ex (x,y), [1, size(x)]));

gradex = @(x,y) cat(1, ...
    reshape (evaluate_fe_function_gradient(exsol.u1,exsol.fespace_u,[x(:)';y(:)']),[1,2,size(x)]), ...
    reshape (evaluate_fe_function_gradient(exsol.u2,exsol.fespace_u,[x(:)';y(:)']),[1,2,size(x)]));

pex = @(x,y) evaluate_fe_function(exsol.p,exsol.fespace_p,[x(:)';y(:)']);

% degree of polynomials for pressure to be used in outflow 
deg = 1;

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

        % compute error
        errsp = sp_l2_error (is{2}.space_p, is{2}.msh, is{2}.p, pex);
        errsu = sp_h1_error (is{2}.space_v, is{2}.msh, is{2}.vel, velex, gradex);
        
        err3(count,freq+1) = errsu + errsp
        
        msh_prc = msh_precompute(is{2}.msh);
        h_iga(count) = max(msh_prc.element_size);
    end
end

save(['data_figure13/errs_P2P1'],'err3','h_iga')

