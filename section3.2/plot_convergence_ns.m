load('data/err')
load('data/h')
n = length(h);

loglog(h(1:end),err(1:end),'.-','Markersize',10)
hold on
expo = 2;
loglog(h,h.^(expo)*min(err(1))/(h(1)^(expo))*2,'--k')
axis([min(h)/1.1 max(h)*1.1 1e-5 1e-2])
set(gca,'fontsize',17)
xticks([0.0078 0.0221 0.0625])
%legend({'Approximate solution','$h^2$'},'Interpreter','latex','Location','Northwest')
xlabel('$h$')
ylabel('$\Vert u - u^h \Vert_{H^1} + \Vert p - p^h \Vert_{L^2}$')
axis square