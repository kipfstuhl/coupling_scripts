close all
clc

load('data_figure13/errs_P2P1.mat')

subplot(1,2,1)

loglog(h_iga,err3,'.-','Markersize',10,'Linewidth',1)
xlim([min(h_iga)/1.1 max(h_iga)*1.1])

load('data_figure13/errs_P3P2.mat')
xlabel('h')
ylabel('error')
subplot(1,2,2)

loglog(h_iga,err3,'.-','Markersize',10,'Linewidth',1)
xlabel('h')
xlim([min(h_iga)/1.1 max(h_iga)*1.1])

legend({'$n_\omega = 0$','$n_\omega = 1$','$n_\omega = 2$','$n_\omega = 3$', ...
        '$n_\omega = 4$','$n_\omega = 5$'},'Interpreter','latex','Location','Southeast')
