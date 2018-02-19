function [mat,rhs,jac,nsys,nus,nps,indices] = build_coupled_system_navier_stokes(fespaces_u,fespaces_p,fun,nu,dirichlet_functions,neumann_functions,domain_connectivity,normals,nbasisfunctions,gausspoints,typebasisfunctions)

nsubdomains = size(fespaces_u,2);

nus = zeros(nsubdomains,1);
nps = zeros(nsubdomains,1);
ns = zeros(nsubdomains,1);

indices = {};

As = {};
rhss = {};

aux2 = 0;
for i = 1:nsubdomains
    nus(i) = size(fespaces_u{i}.nodes,1);
    nps(i) = size(fespaces_p{i}.nodes,1);
    ns(i) = 2*nus(i) + nps(i);
    aux1 = aux2 + 1;
    aux2 = aux2 + ns(i);
    indices{end+1} = aux1:aux2;
    
    [As{end+1},rhss{end+1}] = assembler_steady_navier_stokes(fespaces_u{i},fespaces_p{i},fun,nu,dirichlet_functions,neumann_functions);
end

% handle the interfaces
ninterfaces = size(domain_connectivity,1);

B = [];
Bt = [];

for i = 1:ninterfaces
    interface_connectivity = domain_connectivity(i,:);
    [Bnew,Btnew] = build_coupling_block_navier_stokes(fespaces_u,nus,nps,ns,normals,interface_connectivity,nbasisfunctions(i),gausspoints,typebasisfunctions);
    
    B = [B; Bnew];
    Bt = [Bt Btnew]; 
end

% create block matrix on diagonal recursively and jacobian
A = @(u) [];
jac_block11 = @(u) [];
for i = 1:nsubdomains
    A = @(u) blkdiag(A(u),As{i}(u(indices{i})));
    jac_block11 = @(u) blkdiag(jac_block11(u),build_jac_navier_stokes(As{i},u(indices{i}),fespaces_u{i}));
end

totalagmul = size(B,1);
zerosp = sparse(totalagmul,totalagmul);
mat = @(u) [A(u) Bt; B zerosp];

zeroslg = zeros(totalagmul,1);

rhs = [];
for i = 1:nsubdomains
    rhs = [rhs;rhss{i}];
end
rhs = [rhs;zeroslg];

jac = @(u) [jac_block11(u) Bt; B zerosp];
nsys = sum(ns) + totalagmul;





