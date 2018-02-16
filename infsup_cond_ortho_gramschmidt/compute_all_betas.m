clear all
clc

poly1 = 'P1';
poly2 = 'P1';
norm1 = 'L2';
norm2 = 'L2';

betaL2L2 = compute_beta(poly1,poly2,norm1,norm2);
save('data/betaL2L2.mat','betaL2L2');

poly1 = 'P2';
poly2 = 'P2';
norm1 = 'H1';
norm2 = 'L2';

betaH1L2 = compute_beta(poly1,poly2,norm1,norm2);
save('data/betaH1L2.mat','betaH1L2');

norm1 = 'H1';
norm2 = 'H1';

betaH1H1 = compute_beta(poly1,poly2,norm1,norm2);
save('data/betaH1H1.mat','betaH1H1');


