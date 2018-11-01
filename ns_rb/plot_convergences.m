clear all
close all
clc

load('hs.mat');
load('erroronlysin.mat');

nc = size(totalerr,1);
close all
for i = 1:nc
    loglog(hs,totalerr(i,:),'.-','Linewidth',2,'Markersize',10);
    hold on
end
loglog(hs,hs.^2*0.04,'--k','Linewidth',2)
title('Only sin, constants with double support');
legend('n_\omega = 0', ...
       'n_\omega = 2', ...
       'n_\omega = 4', ...
       'n_\omega = 6', ...
       'n_\omega = 8', ...
       'n_\omega = 10','Location','Northeastoutside');
xlabel('$h$');
ylabel('error');
set(gca,'Fontsize',15)
xlim([min(hs)/1.1 max(hs)*1.1])

figure

load('erroronlysinsinglesupport.mat');

nc = size(totalerr,1);
for i = 1:nc
    loglog(hs,totalerr(i,:),'.-','Linewidth',2,'Markersize',10);
    hold on
end
loglog(hs,hs.^2*0.04,'--k','Linewidth',2)
title('Only sin, constants with single support');
legend('n_\omega = 0', ...
       'n_\omega = 2', ...
       'n_\omega = 4', ...
       'n_\omega = 6', ...
       'n_\omega = 8', ...
       'n_\omega = 10','Location','Northeastoutside');
xlabel('$h$');
ylabel('error');
set(gca,'Fontsize',15)
xlim([min(hs)/1.1 max(hs)*1.1])

figure

load('erroronlysincosdoublesupport.mat');

nc = 6;
for i = 1:nc
    loglog(hs,totalerr(i,:),'.-','Linewidth',2,'Markersize',10);
    hold on
end
loglog(hs,hs.^2*0.04,'--k','Linewidth',2)
title('Sin and cos, constants with double support');
legend('n_\omega = 0', ...
       'n_\omega = 1', ...
       'n_\omega = 2', ...
       'n_\omega = 3', ...
       'n_\omega = 4', ...
       'n_\omega = 5','Location','Northeastoutside');
xlabel('$h$');
ylabel('error');
set(gca,'Fontsize',15)
xlim([min(hs)/1.1 max(hs)*1.1])

figure

load('erroronlysincosdoublesupportlobatto.mat');

nc = 6;
for i = 1:nc
    loglog(hs,totalerr(i,:),'.-','Linewidth',2,'Markersize',10);
    hold on
end
loglog(hs,hs.^2*0.04,'--k','Linewidth',2)
title('Sin and cos, constants with double support,Lobatto');
legend('n_\omega = 0', ...
       'n_\omega = 1', ...
       'n_\omega = 2', ...
       'n_\omega = 3', ...
       'n_\omega = 4', ...
       'n_\omega = 5','Location','Northeastoutside');
xlabel('$h$');
ylabel('error');
set(gca,'Fontsize',15)
xlim([min(hs)/1.1 max(hs)*1.1])

figure

load('erroronlysindoublesupport_lag.mat');

nc = 5;
for i = 1:nc
    semilogy(0:2:10,totalerr_lag(:,i),'.-','Linewidth',2,'Markersize',10);
    hold on
end
title('$\lambda$, sin, constants with double support');
legend(['h = ',num2str(hs(1))], ...
       ['h = ',num2str(hs(2))], ...
       ['h = ',num2str(hs(3))], ...
       ['h = ',num2str(hs(4))], ...
       ['h = ',num2str(hs(5))],'Location','Northeastoutside');
xlabel('$n_\omega$');
ylabel('error');
set(gca,'Fontsize',15)
xlim([0 10])

figure

load('erroronlysindoublesupportlobatto_lag.mat');

nc = 5;
for i = 1:nc
    semilogy(0:2:10,totalerr_lag(:,i),'.-','Linewidth',2,'Markersize',10);
    hold on
end
title('$\lambda$, sin, constants with double support, Lobatto');
legend(['h = ',num2str(hs(1))], ...
       ['h = ',num2str(hs(2))], ...
       ['h = ',num2str(hs(3))], ...
       ['h = ',num2str(hs(4))], ...
       ['h = ',num2str(hs(5))],'Location','Northeastoutside');
xlabel('$n_\omega$');
ylabel('error');
set(gca,'Fontsize',15)
xlim([0 10])

figure

load('erroronlysinsinglesupport_lag.mat');

nc = 5;
for i = 1:nc
    semilogy(0:2:10,totalerr_lag(:,i),'.-','Linewidth',2,'Markersize',10);
    hold on
end
title('$\lambda$, sin, constants with single support');
legend(['h = ',num2str(hs(1))], ...
       ['h = ',num2str(hs(2))], ...
       ['h = ',num2str(hs(3))], ...
       ['h = ',num2str(hs(4))], ...
       ['h = ',num2str(hs(5))],'Location','Northeastoutside');
xlabel('$n_\omega$');
ylabel('error');
set(gca,'Fontsize',15)
xlim([0 10])

figure

load('erroronlysincosdoublesupport_lag.mat');

nc = 5;
for i = 1:nc
    semilogy(0:1:5,totalerr_lag(:,i),'.-','Linewidth',2,'Markersize',10);
    hold on
end
title('$\lambda$, sin and cos, constants with double support');
legend(['h = ',num2str(hs(1))], ...
       ['h = ',num2str(hs(2))], ...
       ['h = ',num2str(hs(3))], ...
       ['h = ',num2str(hs(4))], ...
       ['h = ',num2str(hs(5))],'Location','Northeastoutside');
xlabel('$n_\omega$');
ylabel('error');
set(gca,'Fontsize',15)
xlim([0 5])

figure

load('erroronlysincosdoublesupportlobatto_lag.mat');

nc = 5;
for i = 1:nc
    semilogy(0:1:5,totalerr_lag(:,i),'.-','Linewidth',2,'Markersize',10);
    hold on
end
title('$\lambda$, sin and cos, constants with double support, Lobatto');
legend(['h = ',num2str(hs(1))], ...
       ['h = ',num2str(hs(2))], ...
       ['h = ',num2str(hs(3))], ...
       ['h = ',num2str(hs(4))], ...
       ['h = ',num2str(hs(5))],'Location','Northeastoutside');
xlabel('$n_\omega$');
ylabel('error');
set(gca,'Fontsize',15)
xlim([0 5])




