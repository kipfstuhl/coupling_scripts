clc

% author: Luca Pegolotti on 20/11/2017

% This script plots the convergence order computed by
% compute_order_convergence_P2P2_confmesh.m and reproduces the plots of
% Section 3.

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex')

if (~exist('data','dir'))
    error('ERROR: "data" directory does not exist! Run compute_order_convergence_P2P2_confmesh.m first');
end

% here we load the errors 
load('data/globalerror.mat');
load('data/brokenerror.mat');

% mesh sizes
h = 1./[20 28 40 56 80 114 160];

% select the ngammas we plot (upper limit is 1 + 2 * 15 = 31, must be odd 
% number)
run pick_gammas.m

% discard rows in brokenerror that we do not want to plot
brokenerror = brokenerror(ngamma,:);

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
h1 = loglog(h,h.^2*7,':k','Linewidth',1)
legend_entries{end+1} = '$h^2$';

% set axis
axis([h(end)/1.1 h(1)*1.1 2e-4 1])
axis square
grid on
xlabel('$h$');
ylabel('H$^1$ error');
set(gca, 'fontsize',15); 

% add legend
legend(legend_entries,'interpreter','latex','Location','Northeastoutside')

% add title
title('Convergence for conforming mesh')
