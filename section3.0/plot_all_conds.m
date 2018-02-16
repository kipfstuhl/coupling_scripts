set(0,'defaulttextinterpreter','latex')

poly = {'P2P2'};

subplot(1,2,1)
for k = 1:1

   load(['data/cond',poly{k},'_ortho.mat']);
   n = size(cond,2);
    for i = 1:n
        figure(k)
        semilogy(cond{i}(:,1),cond{i}(:,2),'.-','Linewidth',1)
        hold on
    end
    title('Orthonormal basis functions')
    axis([1 31 1e3 1e20])
end
xlabel('$n_\Gamma$')
ylabel('Condition number')


subplot(1,2,2)
for k = 1:1

   load(['data/cond',poly{k},'_nonortho.mat']);
   n = size(cond,2);
    for i = 1:n
        figure(k)
        semilogy(cond{i}(:,1),cond{i}(:,2),'.-','Linewidth',1)
        hold on
    end
    title('Non-orthonormal basis functions')
    axis([1 31 1e3 1e20])
end
xlabel('$n_\Gamma$')
legend('$h = 1/20$','$h = 1/28$','$h = 1/40$','$h = 1/56$','$h = 1/80$','interpreter','latex', ...
       'Location','Northeastoutside')


