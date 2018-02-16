clear all
close all
clc

co = [0 0 1;
      0 0.5 0;
      1 0 0;
      0 0.75 0.75;
      0.75 0 0.75;
      0.75 0.75 0;
      0.25 0.25 0.25];
set(groot,'defaultAxesColorOrder',co)

set(0,'defaulttextinterpreter','latex')

load('data/betaH1L2_cond.mat');

beta = betaH1L2;

n = size(beta,2);

for i = 1:n
    semilogy(beta{i}(:,1),beta{i}(:,2),'Linewidth',1)
    hold on
end
axis([1 321 0 1])

xlabel('$n_\Gamma$')
ylabel('$\widetilde{\beta}$')
legend({'$h = 1/20$','$h = 1/28$','$h = 1/40$','$h = 1/56$','$h = 1/80$','$h = 1/114$','$h = 1/160$'},'Interpreter','latex','Location','Northeastoutside')
axis square
