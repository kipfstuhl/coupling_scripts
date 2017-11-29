clear all
close all
clc

% author: Luca Pegolotti on 21/11/2017

% This script plots the convergence order computed by
% compute_order_convergence_P2P2_nonconfmesh.m and reproduces the plots of
% Section 3.

% NOTE: the script is actually the same as
% plot_order_convergence_P2P2_nonconfmesh.m exept for little details. For
% this reason, the code is cleared from unnecessary comments.

set(0,'defaulttextinterpreter','latex')

if (~exist('data','dir'))
    error('ERROR: "data" directory does not exist! Run compute_order_convergence_P2P2_confmesh.m first');
end

load('data/globalerror.mat');
load('data/brokenerror_nonconforming.mat');

h = 1./[20 28 40 56 80 114 160];

% select the ngammas we plot (upper limit is 1 + 2 * 15 = 31, must be odd 
% number)
ngamma = [1 3 5 7 13 21 25];

brokenerror_nonconforming = brokenerror_nonconforming((ngamma-1)/2+1,:);

nfreq = size(brokenerror_nonconforming,1);

legend_entries = {};
for i = 1:nfreq
    loglog(h,brokenerror_nonconforming(i,:),'.-','Linewidth',1,'Markersize',10)
    hold on
    legend_entries{end+1} = ['$N_\Gamma$ = ',num2str(ngamma(i))];
end

loglog(h,globalerror,'--k','Linewidth',1)
legend_entries{end+1} = 'Global solution';

loglog(h,h.^2*7,':k','Linewidth',1)
legend_entries{end+1} = '$h^2$';

axis([h(end)/1.1 h(1)*1.1 2e-4 1])
axis square
grid on
xlabel('$h$');
ylabel('H$^1$ error');
set(gca, 'fontsize',15); 

legend(legend_entries,'interpreter','latex','Location','Northeastoutside')

title('Convergence for non-conforming mesh')
