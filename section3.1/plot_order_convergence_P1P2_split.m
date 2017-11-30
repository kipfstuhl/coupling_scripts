clc

% author: Luca Pegolotti on 29/11/2017

% This script plots the convergence order of the global error computed by
% compute_order_convergence_P1P2.m 

% we set the interpreter for strings to latex
set(0,'defaulttextinterpreter','latex')

errortype = 'H1';

if (~exist('data','dir'))
    error('ERROR: "data" directory does not exist! Run compute_order_convergence_P2P2_confmesh.m first');
end

% here we load the errors 
load(['data/globalerror_P1_',errortype,'.mat']);
load(['data/error_P1_',errortype,'.mat']);
load(['data/error_P2_',errortype,'.mat']);

% mesh sizes
h = 1./[20 28 40 56 80 114 160];

% select the ngammas we plot (upper limit is 1 + 2 * 15 = 31, must be odd 
% number)
run pick_gammas.m

% discard rows in brokenerror that we do not want to plot
error1 = error1((ngamma-1)/2+1,:);
error2 = error2((ngamma-1)/2+1,:);

% number of frequencies
nfreq = size(error1,1);

legend_entries = {};

% plot convergence error for domain with P1 basis functions
for i = 1:nfreq
    loglog(h,error1(i,:),'.-','Linewidth',1,'Markersize',10)
    hold on
    legend_entries{end+1} = ['$N_\Gamma$ = ',num2str(ngamma(i))];
end

% plot error of global solution for comparison
loglog(h,globalerror,'--k','Linewidth',1)
legend_entries{end+1} = 'Global solution';

% plot h^2 slope
h1 = loglog(h,h*7,':k','Linewidth',1)
legend_entries{end+1} = '$h^1$';

% plot h^2 slope
h1 = loglog(h,h.^2*7,'-.k','Linewidth',1)
legend_entries{end+1} = '$h^2$';

cfig = gca;
cfig.ColorOrderIndex = 1;

% plot convergence error for domain with P1 basis functions
% for i = 1:nfreq
%     loglog(h,error2(i,:),'.--','Linewidth',1,'Markersize',10)
%     hold on
% end



% set axis
axis([h(end)/1.1 h(1)*1.1 3e-4 1])
axis square
grid on
xlabel('$h$');
ylabel('H$^1$ error');
set(gca, 'fontsize',15); 

% add legend
legend(legend_entries,'interpreter','latex','Location','Northeastoutside')

% add title
title('Convergence for conforming mesh')
