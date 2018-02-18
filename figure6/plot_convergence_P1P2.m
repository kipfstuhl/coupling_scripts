clear all
clc

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex')

N = [20 28 40 56 80 114 160];
h = 1./N;
nrefs = length(N);

gammas = [1 3 5 7 9 13];

leglab = {};

load('data_figure6/globalerrorP1.mat');

subplot(1,2,1)
load('data_figure6/error_domain1.mat');

for i = gammas
    freq = (i + 1)/2;
    loglog(h,error1(freq,:),'.-','Linewidth',1,'Markersize',10);
    leglab{end+1} = ['$n_\Gamma = ',num2str(i),'$'];
    hold on
end

loglog(h,globalerror,'--k','Linewidth',1);

xlabel('$n_\Gamma$');
ylabel('error');
legend(leglab,'Location','Southeast','interpreter','latex');
set(gca,'Fontsize',13)
axis([min(h)/1.1 max(h)*1.1 1e-4 1])
title('$\Omega_1$')

load('data_figure6/globalerrorP2.mat');

subplot(1,2,2)
load('data_figure6/error_domain2.mat');

for i = gammas
    freq = (i + 1)/2;
    loglog(h,error2(freq,:),'.-','Linewidth',1,'Markersize',10);
    leglab{end+1} = ['$n_\Gamma = ',num2str(i),'$'];
    hold on
end

loglog(h,globalerror,'--k','Linewidth',1);

xlabel('$n_\Gamma$');
legend(leglab,'Location','Southeast','interpreter','latex');
set(gca,'Fontsize',13)
axis([min(h)/1.1 max(h)*1.1 1e-4 1])
title('$\Omega_2$')
