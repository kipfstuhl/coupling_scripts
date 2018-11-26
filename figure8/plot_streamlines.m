clear all
close all
clc

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex')

% position of the lines for the partition of the domains
xline1 = 0.15;
yline1 = 0.15;

xline2 = 0.7;
yline2 = 0.35;

c = [102/255, 0/255, 1;
  0 0.5 0;
  0 0.75 0.75;
  1, 0, 0;
  0.75 0 0.75;
  0.75 0.75 0;
  0.25 0.25 0.25];

exsolstr = load('data_figure8/exact_solution.mat');

exsol = exsolstr.sol;

subplot(2,2,1)
streamlines_cavity(exsol,'k',1);

fig = gcf;

% change linestile of exact solution to dashed
lines = fig.Children(1).Children;
for i = 1:length(lines)
    if norm(lines(i).Color - [0 0 0]) == 0
        lines(i).LineStyle = '--';
    end
end

subplot(2,2,2)

% plot single points with colors for legend (they will not be
% displayed)

for i = 1:3
    plot(0,0,'color',c(i,:),'Linewidth',1);
    hold on
end

streamlines_cavity(exsol,'k',1);


lines = fig.Children(1).Children;
for i = 1:length(lines)
    if norm(lines(i).Color - [0 0 0]) == 0
        lines(i).LineStyle = '--';
    end
end

subplot(2,2,3)
streamlines_cavity(exsol,'k',1);

lines = fig.Children(1).Children;
for i = 1:length(lines)
    if norm(lines(i).Color - [0 0 0]) == 0
        lines(i).LineStyle = '--';
    end
end

subplot(2,2,4)
streamlines_cavity(exsol,'k',1);

lines = fig.Children(1).Children;
for i = 1:length(lines)
    if norm(lines(i).Color - [0 0 0]) == 0
        lines(i).LineStyle = '--';
    end
end

count = 0;

for n = [16 32 128]

    count = count + 1;
    load(['data_figure8/intsol_n',num2str(n),'.mat'])
    load(['data_figure8/intsoldom1_n',num2str(n),'.mat'])
    load(['data_figure8/intsoldom2_n',num2str(n),'.mat'])
    load(['data_figure8/intsoldom4_n',num2str(n),'.mat'])
    
    % plot streamlines in the top left
    figure(1)
    subplot(2,2,1)
    streamlines_cavity(intsol,c(count,:),1);
    hold on
    axis([0 0.15 0.85 1])
    box
    axis square

    % plot streamlines in the top right
    subplot(2,2,2)
    streamlines_cavity(intsol,c(count,:),1);
    hold on
    axis([0.85 1 0.85 1])
    box
    axis square

    % plot streamlines in the bottom left
    subplot(2,2,3)
    streamlines_cavity(intsol,c(count,:),1);
    hold on
    axis([0 0.25 0 0.25])
    box
    axis square

    % plot streamlines in the bottom right
    subplot(2,2,4)
    streamlines_cavity(intsol,c(count,:),1);
    hold on
    axis([0.6 1 0 0.4])
    box
    axis square
end

subplot(2,2,3)
plot([0 xline1],[yline1 yline1],'--b')
plot([0 0.3],[yline2 yline2],'--b')
plot([xline1 xline1],[0 0.3],'--b')
plot([xline2 xline2],[0 yline2],'--b')

subplot(2,2,4)
plot([0 xline1],[yline1 yline1],'--b')
plot([0.55 1],[yline2 yline2],'--b')
plot([xline1 xline1],[0 yline2],'--b')
plot([xline2 xline2],[0 yline2],'--b')

subplot(2,2,2)
legend({'$h = 1/16$','$h = 1/32$','$h = 1/128$','Fine solution'},'interpreter','latex');