function [Atot,B_coupling,J,rhs,nurb,int_dofs,vel,Ndofs,space_v,space_p,geometry,msh] =  assemble_system_iga_inflow (method_data,ls,nfreqs)
create_inflow_nurb
geo_name = nurb;

% boundary conditions
drchlt_sides = [1 3 4];

viscosity = @(x, y) ones (size (x));

fx = @(x, y) 0*x;
fy = @(x, y) 0*y;

f  = @(x, y) cat(1, ...
    reshape (fx (x,y), [1, size(x)]), ...
    reshape (fy (x,y), [1, size(x)]));
r_in = 0.5;
h  = @(x,y,iside) cat(1, ...
    reshape(-((y/r_in).^2 - 1)*(iside == 1),[1 size(x)]), ...
    zeros([1 size(x)]));
data_names = fieldnames (method_data);
for iopt  = 1:numel (data_names)
    eval ([data_names{iopt} '= method_data.(data_names{iopt});']);
end

% load geometry
geometry = geo_load (geo_name);

% Compute the mesh structure using the finest mesh
switch (upper(element_name))
    case {'RT', 'TH', 'NDL'}
        [~, zeta] = kntrefine (geometry.nurbs.knots, nsub-1, degree, regularity);
    case {'SG'}
        [~, zeta] = kntrefine (geometry.nurbs.knots, 2*nsub-1, degree, regularity);
end
rule       = msh_gauss_nodes (nquad);
[qn, qw]   = msh_set_quad_nodes (zeta, rule);
msh        = msh_cartesian (zeta, qn, qw, geometry);

% Compute the space structures
[space_v, space_p] = sp_bspline_fluid (element_name, ...
    geometry.nurbs.knots, nsub, degree, regularity, msh);

% Assemble the matrices
if (msh.rdim == 2)
    fun_one = @(x, y) ones (size(x));
elseif (msh.rdim == 3)
    fun_one = @(x, y, z) ones (size(x));
end

restrict = @(A,dim1,dim2) A(dim1,dim2);
extend = @(vel,sol_int,int_dofs) extend_solution(vel,sol_int,int_dofs);

A = op_gradu_gradv_tp (space_v, space_v, msh, viscosity);
B = op_div_v_q_tp (space_v, space_p, msh);
C = @(u) op_ugradu_v_tp (space_v, space_v, msh, u);
E = op_f_v_tp (space_p, msh, fun_one).';
F = op_f_v_tp (space_v, msh, f);

vel   = zeros (space_v.ndof, 1);
press = zeros (space_p.ndof, 1);

nbasis = 2*((2*nfreqs(1)+1) + (2*nfreqs(2)+1) + (2*nfreqs(3)+1));

B_coupling = sparse(space_v.ndof + space_p.ndof,nbasis);

vert_x = 1.75;
vert_y = 0;
x_p = [vert_x vert_y];

count = 1;
interface_side = 2;
for j = 1:3 % interface
    msh_side = msh_eval_boundary_side (msh, interface_side);
    if (strcmpi (element_name, 'RT') || strcmpi (element_name, 'NDL'))
        msh_side_from_interior = msh_boundary_side_from_interior (msh, interface_side);
        
        sp_side = space_v.constructor (msh_side_from_interior);
        sp_side = struct (sp_precompute (sp_side, msh_side_from_interior, 'value', true));
        sp_side.dofs = 1:sp_side.ndof;
    else
        sp_side  = sp_eval_boundary_side (space_v, msh_side);
    end
    
    x = cell (msh_side.rdim, 1);
    for idim = 1:msh_side.rdim
        x{idim} = squeeze (msh_side.geo_map(idim,:,:));
    end
    
    for i = 0:nfreqs(j)
        if (i == 0)
            bx = zeros(space_v.ndof + space_p.ndof,1);
            by = zeros(space_v.ndof + space_p.ndof,1);
            if (j == 1)
                f_b = @(x,y) y < vert_y;
                g = @(x,y,iside) cat(1,reshape(f_b(x,y),[1,size(y)]),zeros([1,size(x)]));
                gval = reshape (g (x{:}, interface_side), msh.rdim, msh_side.nqn, msh_side.nel);
                bx(sp_side.dofs) = op_f_v (sp_side, msh_side, gval);
                
                g = @(x,y,iside) cat(1,zeros([1,size(x)]),reshape(f_b(x,y),[1,size(y)]));
                gval = reshape (g (x{:}, interface_side), msh.rdim, msh_side.nqn, msh_side.nel);
                by(sp_side.dofs) = op_f_v (sp_side, msh_side, gval);
            elseif (j == 3)
                f_b = @(x,y) -(y > vert_y);
                g = @(x,y,iside) cat(1,reshape(f_b(x,y),[1,size(y)]),zeros([1,size(x)]));
                gval = reshape (g (x{:}, interface_side), msh.rdim, msh_side.nqn, msh_side.nel);
                bx(sp_side.dofs) = op_f_v (sp_side, msh_side, gval);
                
                g = @(x,y,iside) cat(1,zeros([1,size(x)]),reshape(f_b(x,y),[1,size(y)]));
                gval = reshape (g (x{:}, interface_side), msh.rdim, msh_side.nqn, msh_side.nel);
                by(sp_side.dofs) = op_f_v (sp_side, msh_side, gval);
            end
        else
            bx = zeros(space_v.ndof + space_p.ndof,2);
            by = zeros(space_v.ndof + space_p.ndof,2);
            if (j == 1)
                f_b = @(x,y,xp,l) sin(i*sqrt((x-xp(1)).^2+(y-xp(2)).^2)/l*pi) .* (y < vert_y);
                
                g = @(x,y,iside) cat(1,reshape(f_b(x,y,x_p,ls(j)),[1,size(y)]),zeros([1,size(x)]));
                gval = reshape (g (x{:}, interface_side), msh.rdim, msh_side.nqn, msh_side.nel);
                bx(sp_side.dofs,1) = op_f_v (sp_side, msh_side, gval);
                
                g = @(x,y,iside) cat(1,zeros([1,size(x)]),reshape(f_b(x,y,x_p,ls(j)),[1,size(y)]));
                gval = reshape (g (x{:}, interface_side), msh.rdim, msh_side.nqn, msh_side.nel);
                by(sp_side.dofs,1) = op_f_v (sp_side, msh_side, gval);
                
                f_b = @(x,y,xp,l) cos(i*sqrt((x-xp(1)).^2+(y-xp(2)).^2)/l*pi) .* (y < vert_y);
                
                g = @(x,y,iside) cat(1,reshape(f_b(x,y,x_p,ls(j)),[1,size(y)]),zeros([1,size(x)]));
                gval = reshape (g (x{:}, interface_side), msh.rdim, msh_side.nqn, msh_side.nel);
                bx(sp_side.dofs,2) = op_f_v (sp_side, msh_side, gval);
                
                g = @(x,y,iside) cat(1,zeros([1,size(x)]),reshape(f_b(x,y,x_p,ls(j)),[1,size(y)]));
                gval = reshape (g (x{:}, interface_side), msh.rdim, msh_side.nqn, msh_side.nel);
                by(sp_side.dofs,2) = op_f_v (sp_side, msh_side, gval);
            elseif (j == 3)
                f_b = @(x,y,xp,l) -sin(i*sqrt((x-xp(1)).^2+(y-xp(2)).^2)/l*pi) .* (y > vert_y);
                
                g = @(x,y,iside) cat(1,reshape(f_b(x,y,x_p,ls(j)),[1,size(y)]),zeros([1,size(x)]));
                gval = reshape (g (x{:}, interface_side), msh.rdim, msh_side.nqn, msh_side.nel);
                bx(sp_side.dofs,1) = op_f_v (sp_side, msh_side, gval);
                
                g = @(x,y,iside) cat(1,zeros([1,size(x)]),reshape(f_b(x,y,x_p,ls(j)),[1,size(y)]));
                gval = reshape (g (x{:}, interface_side), msh.rdim, msh_side.nqn, msh_side.nel);
                by(sp_side.dofs,1) = op_f_v (sp_side, msh_side, gval);
                
                f_b = @(x,y,xp,l) -cos(i*sqrt((x-xp(1)).^2+(y-xp(2)).^2)/l*pi) .* (y > vert_y);
                
                g = @(x,y,iside) cat(1,reshape(f_b(x,y,x_p,ls(j)),[1,size(y)]),zeros([1,size(x)]));
                gval = reshape (g (x{:}, interface_side), msh.rdim, msh_side.nqn, msh_side.nel);
                bx(sp_side.dofs,2) = op_f_v (sp_side, msh_side, gval);
                
                g = @(x,y,iside) cat(1,zeros([1,size(x)]),reshape(f_b(x,y,x_p,ls(j)),[1,size(y)]));
                gval = reshape (g (x{:}, interface_side), msh.rdim, msh_side.nqn, msh_side.nel);
                by(sp_side.dofs,2) = op_f_v (sp_side, msh_side, gval);
            end
        end
        if (i == 0)
            B_coupling(:,count) = bx;
            count = count + 1;
            B_coupling(:,count) = by;
            count = count + 1;
        else
            B_coupling(:,[count count+1]) = bx;
            count = count + 2;
            B_coupling(:,[count count+1]) = by;
            count = count + 2;
        end
    end
end

% Apply Dirichlet  boundary conditions. For RT and NDL elements the normal
%  component is imposed strongly, and the tangential one is imposed weakly.
if (strcmpi (element_name, 'RT') || strcmpi (element_name, 'NDL'))
    [N_mat, N_rhs] = ...
        sp_weak_drchlt_bc (space_v, msh, drchlt_sides, h, viscosity, Cpen);
    [vel_drchlt, drchlt_dofs] = sp_drchlt_l2_proj_udotn (space_v, msh, drchlt_sides, h);
else
    [vel_drchlt, drchlt_dofs] = sp_drchlt_l2_proj (space_v, msh, h, drchlt_sides);
    N_mat = sparse (space_v.ndof, space_v.ndof, 1);
    N_rhs = zeros (space_v.ndof, 1);
end

vel(drchlt_dofs) = vel_drchlt;
int_dofs = setdiff (1:space_v.ndof, drchlt_dofs);
rhs_dir  = @(u) -A(int_dofs, drchlt_dofs)*vel(drchlt_dofs) + N_mat(int_dofs,drchlt_dofs)*vel(drchlt_dofs) - restrict(C(extend(vel,u,int_dofs)),int_dofs,drchlt_dofs)*vel(drchlt_dofs);
Ndofs = numel(int_dofs) + space_p.ndof;

% With natural boundary condition, the constraint on the pressure is not needed.
Atot =   @(u) [ A(int_dofs, int_dofs) - N_mat(int_dofs, int_dofs) + restrict(C(extend(vel,u,int_dofs)),int_dofs,int_dofs), -B(:,int_dofs).';
    -B(:,int_dofs), sparse(size (B,1), size (B,1))];
rhs = @(u) [F(int_dofs) + N_rhs(int_dofs) + rhs_dir(u);
    B(:, drchlt_dofs)*vel(drchlt_dofs)];
C_jac = @(u) op_ugradu_v_jac_tp (space_v, space_v, msh, u);
J =@(u) [ A(int_dofs, int_dofs) - N_mat(int_dofs, int_dofs) + restrict(C_jac(extend(vel,u,int_dofs)),int_dofs,int_dofs), -B(:,int_dofs).';
    -B(:,int_dofs), sparse(size (B,1), size (B,1))];
B_coupling(drchlt_dofs,:) = [];

function vel = extend_solution(vel,sol_int,int_dofs)
vel(int_dofs) = sol_int(1:numel(int_dofs));


