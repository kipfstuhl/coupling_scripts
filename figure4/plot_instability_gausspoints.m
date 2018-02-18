clear all
clc

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex')

N = [20 28 40 56 80 114 160];
h = 1./N;
nrefs = length(N);

gammas = 1:2:31;

leglab = {};

subplot(1,2,1)
load('data_figure4/brokenerror_nonconf_2gausspoints.mat');

for i = 1:nrefs
    semilogy(gammas,brokenerror(:,i),'.-','Linewidth',1,'Markersize',10);
    leglab{end+1} = ['$h = 1/',num2str(N(i)),'$'];
    hold on
end

xlabel('$n_\Gamma$');
legend(leglab,'Location','Northeast','interpreter','latex');
set(gca,'Fontsize',13)
axis([1 31 1e-4 1])
title('2 Gauss quadrature nodes')

subplot(1,2,2)
load('data_figure4/brokenerror_nonconf.mat');

for i = 1:nrefs
    semilogy(gammas,brokenerror(:,i),'.-','Linewidth',1,'Markersize',10);
    leglab{end+1} = ['$h = 1/',num2str(N(i)),'$'];
    hold on
end

xlabel('$n_\Gamma$');
ylabel('error');
legend(leglab,'Location','Northeast','interpreter','latex');
set(gca,'Fontsize',13)
axis([1 31 1e-4 1])
title('4 Gauss quadrature nodes')