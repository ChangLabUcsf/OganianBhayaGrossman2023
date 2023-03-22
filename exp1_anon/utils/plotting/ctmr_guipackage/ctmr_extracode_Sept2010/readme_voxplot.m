
Tthreshold=4.5;% 
Dthreshold=8;% distance threshold from surface in mm

% optional: you may want to only use the voxels from the t-map that fall within the 
% gray matter: reslice the gray matter to the space of the functional scans
% and set a threshold. 0.1 may work
GMthreshold=0.1;

% calculate surface T values
[xyz,t_surf] = voxplot_func_gm(Tthreshold,Dthreshold,GMthreshold);

% render on image
ssize=2;
v='r';
ctmr_vox_plot(xyz,t_surf,ssize,v)
colorbar

%%

Tthreshold=4.5;% 
Dthreshold=8;% distance threshold from surface in mm

% optional: you may want to only use the voxels from the t-map that fall within the 
% gray matter: reslice the gray matter to the space of the functional scans
% and set a threshold. 0.1 may work
GMthreshold=0.1;

% calculate surface T values
[xyz,t_surf] = voxplot_func_gm(Tthreshold,Dthreshold,[]);

%% render on image
ssize=2;
v='r';
ctmr_vox_plot(xyz,t_surf,ssize,v)
colorbar
