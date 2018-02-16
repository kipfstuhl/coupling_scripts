function plot_streamlines_cavity(sol,color,res)

xcenter = 0.546;
ycenter = 0.594;

xs = [0.01 0.03 0.04 0.05];
ys = [0.01 0.03 0.04 0.05];

plot_closed_streamlines(sol,xs,ys,0.05,0.06,color,res);

hold on  

% resolve bottom right eddy

xs = [0.98 0.96 0.94 0.92 0.9 0.9];
ys = [0.02 0.04 0.06 0.08 0.1 0.12];

plot_closed_streamlines(sol,xs,ys,0.88,0.1,color,res);

% resolve central part

xs = [0.06 0.1 0.2 0.3 0.4 0.5 0.55];
ys = xs.^0*0.6;

%plot(xs,ys,'x')

np = length(xs);

plot_closed_streamlines(sol,xs,ys,xcenter,ycenter,color,res);

% resolve outer lines

ys = [0.999];
xs = ys.^0*0.44;

plot_closed_streamlines(sol,xs,ys,xcenter,ycenter,color,res);

ys = [0.99];
xs = [0.995];

plot_closed_streamlines(sol,xs,ys,xcenter,ycenter,color,res);

axis([0 1 0 1])
axis square

