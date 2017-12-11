function blocks = couple_navier_stokes_solutions(fespaces_u,fespaces_p, bcs_flags,base_freq,niterations,label,gausspoints)


ncouplings = size(fespaces_u,2);

blocks = {};
nnodes_u = zeros(ncouplings,1);
nnodes_p = zeros(ncouplings,1);

if (strcmp(label,'xpar'))
    index_com = 1;
elseif (strcmp(label,'ypar'))
    index_com = 2;
else
    error('Label not recognized!');
end

for i = 1:ncouplings
    blocks{end+1} = [];
    nnodes_u(i) = size(fespaces_u{i}.nodes,1);
    nnodes_p(i) = size(fespaces_p{i}.nodes,1);
end

for i = 1:niterations
    bsx = {};
    bsy = {};
    bsp = {};
    
    for j = 1:ncouplings
        bsx{end+1} = zeros(nnodes_u(j),1);
        bsy{end+1} = zeros(nnodes_u(j),1);
        bsp{end+1} = zeros(nnodes_p(j),1);
    end
    
    freq = i - 1;
    
    if (i == 1)
        for j = 1:ncouplings
            bsx{j} = apply_neumann_bc(bsx{j},fespaces_u{j},@(x) bcs_flags(:,j),gausspoints);
            blocks{j} = [blocks{j};[bsx{j};bsy{j};bsp{j}]'];
            blocks{j} = [blocks{j};[bsy{j};bsx{j};bsp{j}]'];
        end
    else
        for j = 1:ncouplings
            bsx{j} = apply_neumann_bc(bsx{j},fespaces_u{j},@(x) bcs_flags(:,j)*sin(x(index_com)*pi*base_freq*freq),gausspoints);
            blocks{j} = [blocks{j};[bsx{j};bsy{j};bsp{j}]'];
            blocks{j} = [blocks{j};[bsy{j};bsx{j};bsp{j}]'];

            bsx{j} = bsx{j}*0;   
         
            bsx{j} = apply_neumann_bc(bsx{j},fespaces_u{j},@(x) bcs_flags(:,j)*cos(x(index_com)*pi*base_freq*freq),gausspoints);
            blocks{j} = [blocks{j};[bsx{j};bsy{j};bsp{j}]'];
            blocks{j} = [blocks{j};[bsy{j};bsx{j};bsp{j}]'];
        end
    end
end