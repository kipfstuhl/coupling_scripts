function [R] = compute_orthonormalization_matrix_onlysin(frequencies,sampling,l,scale)
x = linspace(0,l,sampling)';
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

% sample non-orthonormal basis
V = zeros(sampling,frequencies+1);
V(:,1) = ones(sampling,1);

for i = 2:frequencies+1
   V(:,i) = sin(x*pi*i*scale/l);
end
 
[~,R] = qr(C*V,0);
