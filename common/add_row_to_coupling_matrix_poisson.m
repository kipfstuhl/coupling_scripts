function [B] = add_row_to_coupling_matrix_poisson(B,fespace,boundary_index,freq,varargin)
ngausspoints = 2;
if (nargin > 4)
    ngausspoints = varargin{1};
end

n = size(fespace.nodes,1);
b = zeros(1,n);

bcflags = [0 0 0 0];
bcflags(boundary_index) = 1;

% in this case the interface is parallel to the x axis
if (boundary_index == 1 || boundary_index == 3)
    L = fespace.mesh.L;
    ii = 1; % index along which we integrate
else
    L = fespace.mesh.H;
    ii = 2; % index along which we integrate
end

if (freq == 0)
    % rows of the coupling matrix are constructed similarly to neumann
    % terms (integral on boundaries of the product of the trace of finite
    % element function and a given function)
    b = apply_neumann_bc(b,fespace,@(x) bcflags);
    B = [B;b];
else
    % integrate sin and cos along the interface
    b = apply_neumann_bc(b,fespace,@(x) bcflags * sin(x(ii)*pi*freq/L),ngausspoints);
    B = [B;b];

    b = b*0;
    b = apply_neumann_bc(b,fespace,@(x) bcflags * cos(x(ii)*pi*freq/L),ngausspoints);
    B = [B;b];
end