clc
clear all

xline1 = 0.15;
yline1 = 0.15;

xline2 = 0.7;
yline2 = 0.35;


nex = [16 23 32 45 64];

leg = {};
for n = nex
    load(['data/intsoldom1_n',num2str(n),'.mat']);

    [x,s] = get_stress(sol1,0.35,'Xpar',1e-8);

    h = plot(x,s,'Linewidth',1);
    hold on
    leg{end+1} = ['nx = ',num2str(n)];
end

legend(leg)

exsolstr = load('data/exact_solution.mat');

exsol = exsolstr.sol;

[x,s] = get_stress(exsol,0.35,'Xpar',1e-8);

plot(x,s,'--k','Linewidth',1);