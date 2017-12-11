clc

% author: Luca Pegolotti on 1/12/2017

% This script plots the betah coefficient computed in
% compute_betaj_P2P2_confmesh.m

set(0,'defaulttextinterpreter','latex')

load('data/betah.mat')
ngamma = 1:2:31;

N = [20 28 40 56 80 114 160];

legends = {};

for i = 1:size(betah,2)
    semilogy(ngamma,abs(betah(:,i)))
    hold on
    legends{end+1} = ['$h$ = 1/',num2str(N(i))];
end

legend(legends,'interpreter','latex')
xlabel('$N_\Gamma$');
ylabel('$\beta_h$');
axis([1 31 1e-9 1])
set(gca,'fontsize',15)
