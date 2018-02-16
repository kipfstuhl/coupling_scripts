clc
close all
% author: Luca Pegolotti on 1/12/2017

% This script plots the betah coefficient computed in
% compute_betaj_P2P2_confmesh.m

set(0,'defaulttextinterpreter','latex')

load('data/betah_h1l2.mat')


for i = 1:size(betah_h1l2,2)
    bh1l2 = betah_h1l2{i};
    n = length(bh1l2);
    bf = 1:2:2*n-1;
    semilogy(bf,bh1l2,'.-','Linewidth',1,'Markersize',10)
    hold on
end

legends = {'h = 1/20', 'h = 1/28', 'h = 1/40', 'h = 1/56', 'h = 1/80', 'h = 1/114', 'h = 1/160'};



figure(1)
legend(legends,'interpreter','latex')
xlabel('$N_\Gamma$');
ylabel('$\beta$');
%axis([1 200 0 1e10])
set(gca,'fontsize',15)
title('H1-L2 norms')
