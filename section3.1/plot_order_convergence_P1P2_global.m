clc

% author: Luca Pegolotti on 29/11/2017

% This script plots the convergence order of the global error computed by
% compute_order_convergence_P1P2.m 

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex')

if (~exist('data','dir'))
    error('ERROR: "data" directory does not exist! Run compute_order_convergence_P2P2_confmesh.m first');
end

% here we load the errors 
load('data/globalerror_P1_H1.mat');
load('data/brokenerror_P2P1_H1.mat');

% mesh sizes
h = 1./[20 28 40 56 80 114 160];

% select the ngammas we plot (upper limit is 1 + 2 * 15 = 31, must be odd 
% number)
run pick_gammas.m

% discard rows in brokenerror that we do not want to plot
brokenerror = brokenerror((ngamma-1)/2+1,:);

% number of frequencies
nfreq = size(brokenerror,1);

legend_entries = {};
% plot the convergence of the error for each frequency w (meaning, each line
% represents the convergence of the error obtained by adding Fourier basis
% functions with frequency lower than w)
for i = 1:nfreq
    loglog(h,brokenerror(i,:),'.-','Linewidth',1,'Markersize',10)
    hold on
    legend_entries{end+1} = ['$N_\Gamma$ = ',num2str(ngamma(i))];
end

% plot error of global solution for comparison
loglog(h,globalerror,'--k','Linewidth',1)
legend_entries{end+1} = 'Global solution';

% plot h^2 slope
h1 = loglog(h,h*7,':k','Linewidth',1)
legend_entries{end+1} = '$h^1$';

% set axis
axis([h(end)/1.1 h(1)*1.1 3e-2 1])
axis square
grid on
xlabel('$h$');
ylabel('H$^1$ error');
set(gca, 'fontsize',15); 

% add legend
legend(legend_entries,'interpreter','latex','Location','Northeastoutside')

% add title
title('Convergence for conforming mesh')
