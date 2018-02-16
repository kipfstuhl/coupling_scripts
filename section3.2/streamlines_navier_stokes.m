clc

c = [150, 150, 150;
     100, 100, 100;
     50, 50, 50;
     0, 0, 0]/255;
 
c = [102/255, 0/255, 1;
  0 0.5 0;
  0 0.75 0.75;
  1, 0, 0;
  0.75 0 0.75;
  0.75 0.75 0;
  0.25 0.25 0.25];
 

x0 = [0.05;0.06];
x1 = [0.6;0.6];
x2 = [0.88;0.12];

xline1 = 0.15;
yline1 = 0.15;

xline2 = 0.7;
yline2 = 0.35;

% nex = [16 32 64 128];
nex = [16 32 64 128];

exsolstr = load('data/exact_solution.mat');

exsol = exsolstr.sol;

x0ex = find_eddy_center(exsol,x0);
x1ex = find_eddy_center(exsol,x1);
x2ex = find_eddy_center(exsol,x2);

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

count = 0;

figure(2)
subplot(1,4,1)
plot_streamlines_cavity(exsol,'k',1);

fig = gcf;

lines = fig.Children(1).Children;
for i = 1:length(lines)
    if norm(lines(i).Color - [0 0 0]) == 0
        lines(i).LineStyle = '--';
    end
end

subplot(1,4,2)
plot_streamlines_cavity(exsol,'k',1);


lines = fig.Children(1).Children;
for i = 1:length(lines)
    if norm(lines(i).Color - [0 0 0]) == 0
        lines(i).LineStyle = '--';
    end
end

subplot(1,4,3)
plot_streamlines_cavity(exsol,'k',1);

lines = fig.Children(1).Children;
for i = 1:length(lines)
    if norm(lines(i).Color - [0 0 0]) == 0
        lines(i).LineStyle = '--';
    end
end

subplot(1,4,4)


plot([0.55 0.56],[0.1 0.1],'Color',c(1,:));
hold on
plot([0.55 0.56],[0.1 0.1],'Color',c(2,:));
plot([0.55 0.56],[0.1 0.1],'Color',c(4,:));
plot([0.55 0.56],[0.1 0.1],'--k','Linewidth',1);

plot_streamlines_cavity(exsol,'k',1);

lines = fig.Children(1).Children;
for i = 1:length(lines)
    if norm(lines(i).Color - [0 0 0]) == 0
        lines(i).LineStyle = '--';
    end
end

count_subplot = 1;

for n = nex

    count = count + 1;
    load(['data/intsol_n',num2str(n),'.mat'])
    load(['data/intsoldom1_n',num2str(n),'.mat'])
    load(['data/intsoldom2_n',num2str(n),'.mat'])
    load(['data/intsoldom4_n',num2str(n),'.mat'])

    load(['data/nsys',num2str(n),'.mat'])
    
    if n == 16 || n == 32
        figure(1)
        subplot(1,2,count_subplot)
        count_subplot = count_subplot + 1;
        plot_streamlines_cavity(exsol,'k',1);
        hold on
    %     plot(x0ex(1),x0ex(2),'.k','Markersize',10)
    %     plot(x1ex(1),x1ex(2),'.k','Markersize',10)
    %     plot(x2ex(1),x2ex(2),'.k','Markersize',10)
    %     

        plot_streamlines_cavity(intsol,'r',1);
        title(['$h = 1/',num2str(n),'$']);
        xlabel('$x$');
        ylabel('$y$');
        box

        plot([0 xline1],[yline1 yline1],'--b')
        plot([0 1],[yline2 yline2],'--b')
        plot([xline1 xline1],[0 yline2],'--b')
        plot([xline2 xline2],[0 yline2],'--b')

        set(gca,'Fontsize',17)
    end

     if n~=64
        % plot streamlines in the top left
        figure(2)
        subplot(1,4,1)
        plot_streamlines_cavity(intsol,c(count,:),1);
        hold on
        axis([0 0.15 0.85 1])
        box
        axis square

        % plot streamlines in the top right
        subplot(1,4,2)
        plot_streamlines_cavity(intsol,c(count,:),1);
        hold on
        axis([0.85 1 0.85 1])
        box
        axis square

        % plot streamlines in the bottom left
        subplot(1,4,3)
        plot_streamlines_cavity(intsol,c(count,:),1);
        hold on
        axis([0 0.25 0 0.25])
        box
        axis square

        % plot streamlines in the bottom right
        subplot(1,4,4)
        plot_streamlines_cavity(intsol,c(count,:),1);
        hold on
        axis([0.6 1 0 0.4])
        box
        axis square
    end
end


% simplify data
figure(1)
fig = gcf;

tol = 0.001;

lines = fig.Children(1).Children;
for i = 1:length(lines)
    xs = lines(i).XData;
    ys = lines(i).YData;
    [xs1, ys1, cerr, tol] = reducem(xs', ys',tol);
    diffsize = length(xs) - length(xs1);
    display(['Diff size = ',num2str(diffsize)])
    lines(i).XData = xs1;
    lines(i).YData = ys1;
end

lines = fig.Children(2).Children;
for i = 1:length(lines)
    xs = lines(i).XData;
    ys = lines(i).YData;
    [xs1, ys1, cerr, tol] = reducem(xs', ys',tol);
    diffsize = length(xs) - length(xs1);
    display(['Diff size = ',num2str(diffsize)])
    lines(i).XData = xs1;
    lines(i).YData = ys1;
end

figure(2)
fig = gcf;

lines = fig.Children(1).Children;
for i = 1:length(lines)
    xs = lines(i).XData;
    ys = lines(i).YData;
    size(xs)
    size(ys)
    [xs, ys, cerr, tol] = reducem(xs', ys',tol);
    
    i1 = find(xs < 0.55);
    xs(i1) = [];
    ys(i1) = [];
    
    i1 = find(ys > 0.45);
    xs(i1) = [];
    ys(i1) = [];

    lines(i).XData = xs;
    lines(i).YData = ys;
    
    if (norm(lines(i).Color - [0 0 0]) ~= 0)
        lines(i).LineWidth = 0.5;
    end
end

lines = fig.Children(2).Children;
for i = 1:length(lines)
    xs = lines(i).XData;
    ys = lines(i).YData;
    
    [xs, ys, cerr, tol] = reducem(xs', ys',tol);
    
    i1 = find(xs > 0.3);
    xs(i1) = [];
    ys(i1) = [];
    
    i1 = find(ys > 0.3);
    xs(i1) = [];
    ys(i1) = [];
    
    lines(i).XData = xs;
    lines(i).YData = ys;
    if (norm(lines(i).Color - [0 0 0]) ~= 0)
        lines(i).LineWidth = 0.5;
    end
end

lines = fig.Children(3).Children;
for i = 1:length(lines)
    xs = lines(i).XData;
    ys = lines(i).YData;
    
    
    i1 = find(xs < 0.6);
    xs(i1) = [];
    ys(i1) = [];
    
    i1 = find(ys < 0.6);
    xs(i1) = [];
    ys(i1) = [];
    
    if (length(xs)~= 0 && length(ys)~=0)
        [xs, ys, cerr, tol] = reducem(xs', ys',tol);
    end
    
    lines(i).XData = xs;
    lines(i).YData = ys;
    if (norm(lines(i).Color - [0 0 0]) ~= 0)
        lines(i).LineWidth = 0.5;
    end
end

lines = fig.Children(4).Children;
for i = 1:length(lines)
    xs = lines(i).XData;
    ys = lines(i).YData;
    
    %[xs, ys, cerr, tol] = reducem(xs', ys',tol);
    
    i1 = find(xs > 0.2);
    xs(i1) = [];
    ys(i1) = [];
    
    i1 = find(ys < 0.8);
    xs(i1) = [];
    ys(i1) = [];
    
    if (length(xs)~= 0 && length(ys)~=0)
        [xs, ys, cerr, tol] = reducem(xs', ys',tol);
    end
    
    lines(i).XData = xs;
    lines(i).YData = ys;
    if (norm(lines(i).Color - [0 0 0]) ~= 0)
        lines(i).LineWidth = 0.5;
    end
end

figure(2)
subplot(1,4,3)
plot([0 xline1],[yline1 yline1],'--b')
plot([0 0.3],[yline2 yline2],'--b')
plot([xline1 xline1],[0 0.3],'--b')
plot([xline2 xline2],[0 yline2],'--b')

subplot(1,4,4)
plot([0 xline1],[yline1 yline1],'--b')
plot([0.55 1],[yline2 yline2],'--b')
plot([xline1 xline1],[0 yline2],'--b')
plot([xline2 xline2],[0 yline2],'--b')

box
legend('$h$ = 1/16', ...
       '$h$ = 1/32', ...
       '$h$ = 1/128', ...
       'Fine solution', 'Location', 'NorthEastOutside')
