function [R] = compute_orthonormalization_matrix(frequencies,sampling)
x = linspace(0,1,sampling)';
h = x(2)-x(1);

% compute mass matrix induced by the sampling on (0,1)
d1 = ones(sampling,1) * 2/3*h;
d1(1) = d1(1)/2;
d1(end) = d1(end)/2;
d2 = ones(sampling,1) * 1/6*h;
%M = spdiags(d) + spdiags(d2,1) + spdiags(d2,-1);
M = spdiags([d2 d1 d2],-1:1,sampling,sampling);
% compute cholesky decomposition of the mass matrix
C = chol(M);

bfs = {};
% sample non-orthonormal basis
V = zeros(sampling,frequencies);
V(:,1) = ones(sampling,1);

bfs{end+1} = @(x) x.^0;
for i = 1:Nfrequencies
   index1 = 2*i;
   index2 = 2*i+1;
    
   V(:,index1) = sin(x*pi*i);
   bfs{end+1} = @(x) sin(x*pi*i);
   V(:,index2) = cos(x*pi*i);
   bfs{end+1} = @(x) cos(x*pi*i);
end

[~,R] = qr(C*V,0);