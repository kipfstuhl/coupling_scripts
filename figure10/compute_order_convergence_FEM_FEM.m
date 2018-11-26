clear all
close all
clc

totalerr = [];

for freqq = 0:5
    disp(['Freq = ',num2str(freqq)]);
    
    % arrays for the error and mesh size
    err = [];
    hs = [];
    
    for ref = 1:6
        
        [sols,h] = solve_system_bifurcation_FEM_FEM(ref,freqq);
        
        hs = [hs;h];

        % compute errors
        
        errsu = zeros(3,1);
        errsp = zeros(3,1);
        
        load(['data_figure10/u1s_ref',num2str(ref),'.mat']);
        load(['data_figure10/u2s_ref',num2str(ref),'.mat']);
        load(['data_figure10/ps_ref',num2str(ref),'.mat']);
        
        for i = 1:3 
            dif1 = sols{i}.u1 - u1s{i};
            e1 = compute_H1_error(sols{i}.fespace_u,dif1,@(x) 0,@(x) [0;0]);
            
            dif2 = sols{i}.u2 - u2s{i};
            e2 = compute_H1_error(sols{i}.fespace_u,dif2,@(x) 0,@(x) [0;0]);
            
            errsu(i) = sqrt(e1^2 + e2^2);
            
            difp = sols{i}.p - ps{i};
            e = compute_L2_error(sols{i}.fespace_p,difp,@(x) 0);
            errsp(i) = e;
        end
        erru = sqrt(errsu'*errsu);
        errp = sqrt(errsp'*errsp);
        err = [err erru + errp]
    end
    totalerr = [totalerr;err]
end

save('data_figure10/total_error.mat','totalerr','hs')