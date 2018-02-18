function [B,Bt] = build_coupling_block(fespaces_u,nus,nps,ns,normals,interface_connectivity,nbasisfunctions,gausspoints,typebasisfunctions)

ndomains = length(ns);

bsx = {};
bsy = {};
bsp = {};

B = [];
Bt = [];

% find the lenght of the interface
mincoor = Inf;
maxcoor = -Inf;
if (strcmp(typebasisfunctions,'fourier'))
    if (sum(interface_connectivity == 1) >= 1)
        indices = find(interface_connectivity == 1);
        for i = 1:indices
            mincoor = min(fespaces_u{i}.mesh.xp,mincoor);
            maxcoor = max(fespaces_u{i}.mesh.xp+fespaces_u{i}.mesh.L,maxcoor);
        end
    elseif (sum(interface_connectivity == 2) >= 1)
        indices = find(interface_connectivity == 2);
        for i = 1:indices
            mincoor = min(fespaces_u{i}.mesh.yp,mincoor);
            maxcoor = max(fespaces_u{i}.mesh.yp+fespaces_u{i}.mesh.H,maxcoor);
        end
    elseif (sum(interface_connectivity == 3) >= 1)
        indices = find(interface_connectivity == 3);
        for i = 1:indices
            mincoor = min(fespaces_u{i}.mesh.xp,mincoor);
            maxcoor = max(fespaces_u{i}.mesh.xp+fespaces_u{i}.mesh.L,maxcoor);
        end
    elseif (sum(interface_connectivity == 4) >= 1)
        indices = find(interface_connectivity == 4);
        for i = 1:indices
            mincoor = min(fespaces_u{i}.mesh.yp,mincoor);
            maxcoor = max(fespaces_u{i}.mesh.yp+fespaces_u{i}.mesh.H,maxcoor);
        end
    end
end
lengthinterface = maxcoor - mincoor;

for i = 1:nbasisfunctions
    block_x = [];
    block_y = [];
    blockt_x = [];
    blockt_y = [];
    
    % these are actually only needed when the basis is of fourier type
    block_xx = [];
    block_yy = [];
    blockt_xx = [];
    blockt_yy = [];
    
    for j = 1:ndomains
        bsx = zeros(nus(j),1);
        bsy = zeros(nus(j),1);
        bsp = zeros(nps(j),1);
        
        if (interface_connectivity(j) == 0)
            block_x = [block_x [bsx;bsy;bsp]'];
            block_y = [block_y [bsx;bsy;bsp]'];
            blockt_x = [blockt_x; [bsx;bsy;bsp]];
            blockt_y = [blockt_y; [bsx;bsy;bsp]];
            
            if (i~=1 && strcmp(typebasisfunctions,'fourier'))
                block_xx = [block_xx [bsx;bsy;bsp]'];
                block_yy = [block_yy [bsx;bsy;bsp]'];
                blockt_xx = [blockt_xx; [bsx;bsy;bsp]];
                blockt_yy = [blockt_yy; [bsx;bsy;bsp]];
            end
            
        else
            bc_flags = [0 0 0 0]';
            bc_flags(interface_connectivity(j)) = normals(j,interface_connectivity(j));
            
            if (mod(interface_connectivity(j),2) == 1)
                coor = 1;
            else
                coor = 2;
            end

            if (strcmp(typebasisfunctions,'polynomial'))  
                bf = @(x) x(coor,:).^(i-1);
                bsx = apply_neumann_bc(bsx,fespaces_u{j},@(x) bf(x)*bc_flags,gausspoints);             
                block_x = [block_x [bsx;bsy;bsp]'];
                block_y = [block_y [bsy;bsx;bsp]'];
                
                bsx = apply_dirichlet_bc_rhs(bsx,fespaces_u{j},@(x) [0;0;0;0]);
                blockt_x = [blockt_x;[bsx;bsy;bsp]];
                blockt_y = [blockt_y;[bsy;bsx;bsp]];
            elseif (strcmp(typebasisfunctions,'fourier'))
                freq = (i-1)/lengthinterface;
                bf = @(x) cos(x(coor,:)*freq*pi);
                bsx = apply_neumann_bc(bsx,fespaces_u{j},@(x) bf(x)*bc_flags,gausspoints);             
                block_x = [block_x [bsx;bsy;bsp]'];
                block_y = [block_y [bsy;bsx;bsp]'];
                
                bsx = apply_dirichlet_bc_rhs(bsx,fespaces_u{j},@(x) [0;0;0;0]);
                blockt_x = [blockt_x;[bsx;bsy;bsp]];
                blockt_y = [blockt_y;[bsy;bsx;bsp]];
                
                if (freq > 0)
                    bsx = bsx*0;
                    
                    bf = @(x) sin(x(coor,:)*freq*pi);
                    bsx = apply_neumann_bc(bsx,fespaces_u{j},@(x) bf(x)*bc_flags,gausspoints);             
                    block_xx = [block_xx [bsx;bsy;bsp]'];
                    block_yy = [block_yy [bsy;bsx;bsp]'];

                    bsx = apply_dirichlet_bc_rhs(bsx,fespaces_u{j},@(x) [0;0;0;0]);
                    blockt_xx = [blockt_xx;[bsx;bsy;bsp]];
                    blockt_yy = [blockt_yy;[bsy;bsx;bsp]];
                end
            else
               error('Type of polynomial basis is not supported'); 
            end
        end
    end
    B = [B; block_x; block_xx; block_y; block_yy];
    Bt = [Bt blockt_x blockt_xx blockt_y blockt_yy];
end

B = sparse(B);
Bt = sparse(Bt);