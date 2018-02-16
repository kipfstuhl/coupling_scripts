clear all
close all
clc

norms = {'L2L2','H1L2','H1H1'};

for k = 1:3

    load(['data/beta',norms{k},'.mat']);
   if (k == 1)
       beta = betaL2L2;
   elseif (k == 2)
       beta = betaH1L2;
   elseif (k == 3)
       beta = betaH1H1;
   end
    n = size(beta,2);

    for i = 1:n
        figure(k)
        semilogy(beta{i}(:,1),beta{i}(:,2))
        hold on
    end
    title([norms{k},' norms']);
end
