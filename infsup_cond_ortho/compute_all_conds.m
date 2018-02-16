clear all
clc

poly1 = 'P1';
poly2 = 'P1';

cond = compute_conds(poly1,poly2);
save('data/condP1P1.mat','cond');

poly1 = 'P1';
poly2 = 'P2';

cond = compute_conds(poly1,poly2);
save('data/condP1P2.mat','cond');

poly1 = 'P2';
poly2 = 'P2';

cond = compute_conds(poly1,poly2);
save('data/condP2P2.mat','cond');