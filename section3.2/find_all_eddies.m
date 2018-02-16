clc

x0 = [0.05;0.06];
x1 = [0.6;0.6];
x2 = [0.88;0.12];

exsolstr = load('data/exact_solution.mat');

exsol = exsolstr.sol;

x0ex = find_eddy_center(exsol,x0);
x1ex = find_eddy_center(exsol,x1);
x2ex = find_eddy_center(exsol,x2);

figure(1)
c = [1 1 255;
     234 16 20;
     154 16 234;
     234 172 16;
     16 234 205]/255

for i = 1:5
   if (n == 1 || n == 3)
        plot(0,0,'color',c(i,:),'Linewidth',1)
   end
   hold on
end

nex = [16 32 64 128];

syss = [];
errs0 = [];
errs1 = [];
errs2 = [];

syss_c = [];
errs0_c = [];
errs1_c = [];
errs2_c = [];

count = 0;

leg = {};
for n = nex
    count = count + 1;
    load(['data/intsol_n',num2str(n),'.mat'])
    load(['data/intsoldom1_n',num2str(n),'.mat'])
    load(['data/intsoldom2_n',num2str(n),'.mat'])
    load(['data/intsoldom4_n',num2str(n),'.mat'])

    load(['data/nsys',num2str(n),'.mat'])
    
    syss = [syss nsys];
    x0app = find_eddy_center(sol4,x0);
    x1app = find_eddy_center(sol1,x1);
    x2app = find_eddy_center(sol2,x2);
    
    errs0 = [errs0 norm(x0app-x0ex)];
    errs1 = [errs1 norm(x1app-x1ex)];
    errs2 = [errs2 norm(x2app-x2ex)];
    
    if (n == 16 || n == 64)
        plot_streamlines_cavity(intsol,c(count,:),1);
        leg{end+1} = ['n = ',num2str(n)];
    end
    load(['data/nsys',num2str(n),'_coarse.mat'])
    load(['data/coarse_sol',num2str(n),'.mat'])

    syss_c = [syss_c nsys];
    x0app = find_eddy_center(sol,x0);
    x1app = find_eddy_center(sol,x1);
    x2app = find_eddy_center(sol,x2);
    
    errs0_c = [errs0_c norm(x0app-x0ex)];
    errs1_c = [errs1_c norm(x1app-x1ex)];
    errs2_c = [errs2_c norm(x2app-x2ex)];
    
end

plot_streamlines_cavity(exsol,'--k',1);
hold on
plot(x0ex(1),x0ex(2),'x','Markersize',10)
plot(x1ex(1),x1ex(2),'x','Markersize',10)
plot(x2ex(1),x2ex(2),'x','Markersize',10)

legend(leg);
figure(2)

loglog(syss_c,errs0,'r');
hold on
loglog(syss_c,errs1,'k');
loglog(syss_c,errs2,'b');

loglog(syss_c,errs0_c,'--r');
hold on
loglog(syss_c,errs1_c,'--k');
loglog(syss_c,errs2_c,'--b');
axis square
set(gca,'Fontsize',15)
xlabel('error')
ylabel('dofs')
