clear all
close all
clc

poly = {'P1P1','P1P2','P2P2'};

for k = 1:3

   load(['data/cond',poly{k},'.mat']);
   n = size(cond,2);
    for i = 1:n
        figure(k)
        semilogy(cond{i}(:,1),cond{i}(:,2))
        hold on
    end
    title([poly{k},' polynomials']);
    ylim([1e3 1e20])
end
