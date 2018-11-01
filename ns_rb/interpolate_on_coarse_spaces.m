clear all
clc

load('solutions/exsol.mat');
exsol = sol;

for i = 1:6
    disp(['Interpolating on refinement ',num2str(i)]);
    mesh_in = read_mesh(['../meshes/refinement',num2str(i),'/inflow_distorted.msh']);
    mesh_out1 = read_mesh(['../meshes/refinement',num2str(i),'/outflow1_distorted.msh']);
    mesh_out2 = read_mesh(['../meshes/refinement',num2str(i),'/outflow2_distorted.msh']);
    
    fespace_us = cell(3,1);
    fespace_ps = cell(3,1);
    
    % fespaces inlet
    fespace_us{1} = create_fespace(mesh_in,'P2',[1 0 0 1 1]);
    fespace_ps{1} = create_fespace(mesh_in,'P1',[]);

    % fespaces outlet1
    fespace_us{2} = create_fespace(mesh_out1,'P2',[1 0 1 0 0]);
    fespace_ps{2} = create_fespace(mesh_out1,'P1',[]);

    % fespaces outlet2
    fespace_us{3} = create_fespace(mesh_out2,'P2',[0 1 0 1 0]);
    fespace_ps{3} = create_fespace(mesh_out2,'P1',[]);
    
    u1s = cell(3,1);
    u2s = cell(3,1);
    ps = cell(3,1);
    
    for j = 1:3
       u1s{j} = interp_on_fespace(exsol.fespace_u,exsol.u1,fespace_us{j});
       u2s{j} = interp_on_fespace(exsol.fespace_u,exsol.u2,fespace_us{j});
       ps{j} = interp_on_fespace(exsol.fespace_p,exsol.p,fespace_ps{j});
    end
    save(['solutions/u1s_ref',num2str(i),'.mat'],'u1s');
    save(['solutions/u2s_ref',num2str(i),'.mat'],'u2s');
    save(['solutions/ps_ref',num2str(i),'.mat'],'ps');
end