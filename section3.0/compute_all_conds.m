clear all
clc

poly1 = 'P2';
poly2 = 'P2';

cond = compute_conds_ortho(poly1,poly2);
save('data/condP2P2_ortho.mat','cond');

cond = compute_conds_nonortho(poly1,poly2);
save('data/condP2P2_nonortho.mat','cond');