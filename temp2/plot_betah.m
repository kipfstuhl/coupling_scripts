clc
close all
% author: Luca Pegolotti on 1/12/2017

% This script plots the betah coefficient computed in
% compute_betaj_P2P2_confmesh.m

set(0,'defaulttextinterpreter','latex')

load('data/betah_l2l2.mat')
load('data/nbfsh.mat')
load('data/betah_h1h1.mat')
load('data/condh.mat')


N = [20]
legends = {};

for i = 1:size(betah_h1h1,2)
    try     
        bbl2l2 = betah_l2l2{i};
        l2size = size(bbl2l2);
        figure(1)
        semilogy(bf(1:l2size),bbl2l2,'.-','Linewidth',1,'Markersize',10)
    catch
    end
    hold on
    bf = nbfsh{i};
    cs = condh{i};
    is = bf(cs<1e99);
    cs = cs(cs<1e99);
    bbh1h1 = betah_h1h1{i};
    figure(2)
    semilogy(is,cs,'.-','Linewidth',1,'Markersize',10)
    hold on
    figure(3)
    semilogy(bf,bbh1h1,'.-','Linewidth',1,'Markersize',10)
    hold on
    %legends{end+1} = ['$h$ = 1/',num2str(N(i))];
end

legends = {'h = 1/20', 'h = 1/28', 'h = 1/40', 'h = 1/56', 'h = 1/80', 'h = 1/114', 'h = 1/160'};



figure(1)
legend(legends,'interpreter','latex')
xlabel('$N_\Gamma$');
ylabel('$\beta$');
%axis([1 200 0 1e10])
set(gca,'fontsize',15)
title('L2-L2 norms')

figure(2)
legend(legends,'interpreter','latex')
xlabel('$N_\Gamma$');
ylabel('$\mathcal K(A)$');
axis([1 30 0 1e20])
set(gca,'fontsize',15)
title('Condition number')

figure(3)
legend(legends,'interpreter','latex')
xlabel('$N_\Gamma$');
ylabel('$\beta$');
%axis([1 200 0 1e10])
set(gca,'fontsize',15)
title('H1-H1 norms')
