clear all
clc

poly1 = 'P2';
poly2 = 'P2';
norm1 = 'H1';
norm2 = 'L2';

betaH1L2 = compute_beta(poly1,poly2,norm1,norm2);
save('data/betaH1L2.mat','betaH1L2');

