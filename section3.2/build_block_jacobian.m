function [J] = build_block_jacobian(As,B,u,fespaces_u,fespaces_p)

nfespaces = size(fespaces_u,2);

us = {};
nsu = zeros(size(nfespaces_u.nodes,1),1);
nsp = zeros(size(nfespaces_p.nodes,1),1);

aux = 0;
for i = 1:nfespaces
    nsu(i) = size(fespaces_u{i}.nodes,1);
    nsp(i) = size(fespaces_p{i}.nodes,1);
    us(i) = u(aux+1:aux+nsu(i)*2);
    aux = aux + 2*nsu(i)+nsp(i);
end

Js = {};

aux = 0;
for i = 1:nfespaces
    currindices = aux+1:aux+nsu(i)*2;
    Js{end+1} = build_jac_navier_stokes(As{i},u(currindices),fespaces_u{i});
end

