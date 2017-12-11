clc

% author: Luca Pegolotti on 30/11/2017

% This script plots the convergence wrt the number of basis functions at the interface
% by using the data computed in compute_order_convergence_P2P2_nonconfmesh.m 

set(0,'defaulttextinterpreter','latex')

if (~exist('data','dir'))
    error('ERROR: "data" directory does not exist! Run compute_order_convergence_P2P2_confmesh.m first');
end

load('data/brokenerror_nonconforming_gaussorder3.mat');

N = [20 28 40 56 80 114 160];
gammas = 1:2:31;

legend_entries = {};
for i = 1:length(N)
    semilogy(gammas,brokenerror_nonconforming(:,i),'.-','Linewidth',1,'Markersize',10)
    hold on
    legend_entries{end+1} = ['$h $ = 1/',num2str(N(i))];
end

axis([1 31 2e-4 1])
axis square
grid on
xlabel('$N_\Gamma$');
ylabel('H$^1$ error');
set(gca, 'fontsize',15); 

legend(legend_entries,'interpreter','latex','Location','Northeastoutside')

title('2 Gauss points')
