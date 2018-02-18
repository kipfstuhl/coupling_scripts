clear all
clc

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex');
N = [20 28 40 56 80 114 160];
nrefs = length(N);

leglab = {};

load('data_figure5/beta.mat');

for i = 1:nrefs
    b = beta{i};
    plot(b(:,1),b(:,2),'.-','Linewidth',1,'Markersize',10);
    leglab{end+1} = ['$h = 1/',num2str(N(i)),'$'];
    hold on
end

xlabel('$n_\Gamma$');
ylabel('$\widetilde \beta$');
legend(leglab,'Location','Southwest','interpreter','latex');
set(gca,'Fontsize',13)
axis([1 31 0 1.5])