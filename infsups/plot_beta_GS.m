close all
load('data_infsups/beta_GS.mat');

set(0,'defaulttextinterpreter','latex')

n = size(beta,2);

for i = 1:n
    plot(beta{i}(:,1),beta{i}(:,2),'Linewidth',1)
    hold on
end
axis([1 31 5e-4 2])

legend({'$h = 1/20$','$h = 1/28$','$h = 1/40$','$h = 1/56$','$h = 1/80$','$h = 1/114$','$h = 1/160$'},'interpreter','latex','Location','Southwest')
xlabel('$N_\Gamma$')
ylabel('$\widetilde \beta$');