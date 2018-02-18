clear all
clc

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex')

load('data_figure7/err.mat')
load('data_figure7/h.mat')

loglog(h,err,'.-b','Linewidth',1,'Markersize',10)

hold on
loglog(h,h.^2*min(err(1))/(h(1)^2)*2,'--k')

axis([min(h)/1.1 max(h)*1.1 1e-5 1e-2])

set(gca,'fontsize',17)
xticks([0.0078 0.0221 0.0625])
xlabel('$h$')
ylabel('$\Vert u - u^h \Vert_{H^1} + \Vert p - p^h \Vert_{L^2}$')
axis square

