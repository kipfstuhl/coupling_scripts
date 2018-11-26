clear all
clc

figure

load('data_figure10/total_error.mat');
nc = size(totalerr,1);
for i = 1:nc
    loglog(hs,totalerr(i,:),'.-','Linewidth',2,'Markersize',10);
    hold on
end
loglog(hs,hs.^2*12,'--k','Linewidth',2)
legend('n_\omega = 0', ...
       'n_\omega = 1', ...
       'n_\omega = 2', ...
       'n_\omega = 3', ...
       'n_\omega = 4', ...
       'n_\omega = 5','Location','Southeast');
xlabel('$h$');
ylabel('error');

set(gca,'Fontsize',15)
xlim([min(hs)/1.1 max(hs)*1.1])
ylim([1e-2 5e-1])
