clear all
clc

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex')

N = [20 28 40 56 80 114 160];
h = 1./N;
nrefs = length(N);

gammas = [1 3 5 7 9 13];

leglab = {};

load('data_figure3/globalerror.mat');

subplot(1,2,1)
load('data_figure3/brokenerror_conf.mat');

for i = gammas
    freq = (i + 1)/2;
    loglog(h,brokenerror(freq,:),'.-','Linewidth',1,'Markersize',10);
    leglab{end+1} = ['$n_\Gamma = ',num2str(i),'$'];
    hold on
end

loglog(h,globalerror,'--k','Linewidth',1);

xlabel('$n_\Gamma$');
ylabel('error');
legend(leglab,'Location','Southeast','interpreter','latex');
set(gca,'Fontsize',13)
axis([min(h)/1.1 max(h)*1.1 1e-4 1])
title('Conforming mesh')

subplot(1,2,2)
load('data_figure3/brokenerror_nonconf.mat');

for i = gammas
    freq = (i + 1)/2;
    loglog(h,brokenerror(freq,:),'.-','Linewidth',1,'Markersize',10);
    hold on
end

loglog(h,globalerror,'--k','Linewidth',1);

xlabel('$n_\Gamma$');
legend(leglab,'Location','Southeast','interpreter','latex');
set(gca,'Fontsize',13)
axis([min(h)/1.1 max(h)*1.1 1e-4 1])
title('Non-conforming mesh')
