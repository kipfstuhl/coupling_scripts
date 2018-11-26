close all
clc

load('errs_thP2P1.mat')
load('errs_thP3P2.mat')

subplot(1,2,1)

loglog(h_iga,err2)

xlim([min(h_iga)/1.1 max(h_iga)*1.1])
ylim([1e-4 2])

subplot(1,2,2)

loglog(h_iga,err1)

xlim([min(h_iga)/1.1 max(h_iga)*1.1])
ylim([1e-4 2])

legend('n_\omega = 0','n_\omega = 1','n_\omega = 2','n_\omega = 3','n_\omega = 4','n_\omega = 5','Location','NortheastOutside')
