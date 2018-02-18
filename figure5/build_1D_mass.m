function M = build_1D_mass(x0,x1,h)
nelem  = (x1-x0)/h;

N = nelem*2 + 1;

I1 = 2*h/15; %function 1 times function1
I2 = -h/30; % function 1 times function 3
I3 = h/15; % function 1 times function 2
I4 = 8*h/15; % function 2 times function 2

M = zeros(N,N);
M(1,1) = I1;
M(1,2) = I2;
M(1,nelem+2) = I3;
M(nelem+2,1) = I3;
for i = 2:nelem 
    M(i,i) = 2*I1;
    M(i,i+1) = I2;
    M(i,i-1) = I2;
    M(i,i+nelem+1) = I3;
    M(i,i+nelem) = I3;
    M(i+nelem+1,i) = I3;
    M(i+nelem,i) = I3;
end
M(nelem+1,nelem+1) = I1;
M(nelem+1,nelem) = I2;
M(nelem+1,2*nelem+1) = I3;
M(2*nelem+1,nelem+1) = I3;

for i = nelem+2:N
    M(i,i) = I4;
end

M = sparse(M);
