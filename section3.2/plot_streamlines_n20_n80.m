clc

xline1 = 0.15;
yline1 = 0.15;

xline2 = 0.7;
yline2 = 0.35;

res = 50;

exsolstr = load('data/exact_solution.mat');

exsol = exsolstr.sol;

subplot(1,2,1)
load('data/coarse_sol.mat');

intsol = interpolate_multiple_solutions({sol},exsol.fespace_u,exsol.fespace_p)

plot_streamlines_cavity(exsol,'k',res);
hold on
plot_streamlines_cavity(intsol,'r',res);
title('Complete with N = 20')
box

load('data/ntsol20.mat');
subplot(1,2,2)
plot_streamlines_cavity(exsol,'k',res);
hold on
plot_streamlines_cavity(intsol,'r',res);
plot([0 xline1],[yline1 yline1],'--b')
plot([0 1],[yline2 yline2],'--b')
plot([xline1 xline1],[0 yline2],'--b')
plot([xline2 xline2],[0 yline2],'--b')
title('Splitted with = 15')
box 
