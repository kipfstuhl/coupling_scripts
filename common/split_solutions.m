function [sols,lm] = split_solutions(sol,fespaces_u,fespaces_p,nus,nps,indices)

nsolutions = size(indices,2);

sols = {};
for i = 1:nsolutions
    solvec = sol(indices{i});

    newsol.u1 = solvec(1:nus(i));
    newsol.u2 = solvec(nus(i)+1:2*nus(i));
    newsol.p  = solvec(2*nus(i)+1:2*nus(i)+nps(i));
    newsol.fespace_u = fespaces_u{i};
    newsol.fespace_p = fespaces_p{i};
    
    sols{end+1} = newsol;
end

last_indices = indices{i};
lm = sol(last_indices(2)+1:end);