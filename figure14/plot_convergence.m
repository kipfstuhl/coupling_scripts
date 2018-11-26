close all
clc

load('data_figure14/errs_P3P2.mat')

subplot(1,2,1)

loglog(h_iga,err1,'.-','Markersize',10,'Linewidth',1)

xlim([min(h_iga)/1.1 max(h_iga)*1.1])
ylim([1e-4 2])
xlabel('h')
ylabel('error')

subplot(1,2,2)

loglog(h_iga,err2,'.-','Markersize',10,'Linewidth',1)

xlim([min(h_iga)/1.1 max(h_iga)*1.1])
ylim([1e-4 2])
xlabel('h')

legend('n_\omega = 0','n_\omega = 1','n_\omega = 2','n_\omega = 3','n_\omega = 4','n_\omega = 5','Location','Southeast')
