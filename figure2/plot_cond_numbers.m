clear all
clc

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex')

N = [20 28 40 56 80];
nrefs = length(N);
ngamma = 1:2:31;
leglab = {};

subplot(1,2,1)
load('data_figure2/condnum_nonortho.mat');

for i = 1:nrefs
    semilogy(ngamma,condnum(:,i),'.-','Linewidth',1,'Markersize',10);
    leglab{end+1} = ['$h = 1/',num2str(N(i)),'$'];
    hold on
end

xlabel('$n_\Gamma$');
ylabel('condition number');
legend(leglab,'Location','Southeast','interpreter','latex');
set(gca,'Fontsize',13)
axis([1 31 0 1e20])

subplot(1,2,2)
load('data_figure2/condnum_ortho.mat');

for i = 1:nrefs
    semilogy(ngamma,condnum(:,i),'.-','Linewidth',1,'Markersize',10);
    hold on
end

xlabel('$n_\Gamma$');
legend(leglab);
legend(leglab,'Location','Northwest','interpreter','latex');
set(gca,'Fontsize',13)
axis([1 31 0 1e10])
