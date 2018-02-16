clear all
clc

beta = compute_beta_GS();
save('data_infsups/beta_GS.mat','beta')

beta = compute_beta_default_fourier();
save('data_infsups/beta_default.mat','beta')





