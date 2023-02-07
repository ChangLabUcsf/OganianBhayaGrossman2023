%     Copyright (C) 2009  D. Hermes & K.J. Miller, Dept of Neurology and Neurosurgery, University Medical Center Utrecht
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
%% preprocessing of scans

% coregister + reslice CT to anatomical MR 
%   in SPM5 reference image = MR, source image = CT
% segment MR 


%% startup SPM, spm functions are used in following scripts


%% 1) generate surface to project electrodes to

% get_mask(subject,degree of smoothing,threshold) 
% e.g. get_mask('name',10,0.3) or get_mask('name',16,0.3)

subject='CM';
get_mask(subject,10,0.3);

% surface is saved as name_surface1_10_03.img
% check the surface in spm or mricron 
% it should be a smooth layer on the gray matter surface

%% 2) localize electrodes from ct
ctmr

%% 3) sort unprojected electrodes
sortElectrodes;
% saves as electrodes_locX;

%% 3) combine sorted electrodes in one file

subject='name';
if isequal(subject,'name')
    elec1=load('./data/electrodes_loc1.mat'); % grid 1
    elec2=load('./data/electrodes_loc2.mat'); % grid 2T
    elecmatrix=[elec1.elecmatrix;...
        elec2.elecmatrix];   
    save('./data/electrodes_loc_all.mat','elecmatrix');
end


%% 4) project electrodes 2 surface
% electrodes2surf(subject,localnorm index,do not project electrodes closer than 3 mm to surface)
subject='name';

% electrodes are projected per grid with electrodes2surf.m
% in this example there were 7 grids

% electrodes2surf(
    % 1: subject
    % 2: number of electrodes local norm for projection (0 = whole grid)
    % 3: 0 = project all electrodes, 1 = only project electrodes > 3 mm
    %    from surface, 2 = only map to closest point (for strips)
    % 4: electrode numbers
    % 5: (optional) electrode matrix.mat (if not given, SPM_select popup)
    % 6: (optional) surface.img (if not given, SPM_select popup)
    % 7: (optional) mr.img for same image space with electrode
    %    positions
% automatically saves:
%       a matrix with projected electrode positions 
%       a nifti image with projected electrodes
% saved as electrodes_onsurface_filenumber_inputnr2

% Trial (SN)
[out_els,out_els_ind]=electrodes2surf(subject,4,1,9:24,'electrodes_loc_all.mat','new_mask20_25p5.img','AP_left_grey.nii');


% 1 (1_5) 
[out_els,out_els_ind]=electrodes2surf(subject,5,1,1:32,'./data/electrodes_loc_all.mat','./data/name_surface1_10_03.img','./data/MR.nii');
% 2 (1_4)
[out_els,out_els_ind]=electrodes2surf(subject,4,1,33:48,'./data/electrodes_loc_all.mat','./data/name_surface1_10_03.img','./data/MR.nii');
% 3 (2_4)
[out_els,out_els_ind]=electrodes2surf(subject,4,1,49:64,'./data/electrodes_loc_all.mat','./data/name_surface1_10_03.img','./data/MR.nii');
% 4 (0_2)
[out_els,out_els_ind]=electrodes2surf(subject,0,2,65:72,'./data/electrodes_loc_all.mat','./data/name_surface1_10_03.img','./data/MR.nii');
% 5 (1_2)
[out_els,out_els_ind]=electrodes2surf(subject,0,2,73:80,'./data/electrodes_loc_all.mat','./data/name_surface1_10_03.img','./data/MR.nii');
% 6 (2_2)
[out_els,out_els_ind]=electrodes2surf(subject,0,2,81:88,'./data/electrodes_loc_all.mat','./data/name_surface1_10_03.img','./data/MR.nii');
% 7 (3_4)
[out_els,out_els_ind]=electrodes2surf(subject,4,1,89:96,'./data/electrodes_loc_all.mat','./data/simo_surface1_10_03.img','./data/ANAT_simo131109_5_1-0001.nii');


%% 5) combine electrode files into one 

load('./data/name_electrodesOnsurface1_5.mat'); % grid 1
elecmatrix=out_els;
load('./data/name_electrodesOnsurface1_4.mat'); % grid 2
elecmatrix=[elecmatrix;out_els];
load('./data/name_electrodesOnsurface2_4.mat'); % grid 3
elecmatrix=[elecmatrix;out_els];
load('./data/name_electrodesOnsurface1_0.mat'); % grid 4
elecmatrix=[elecmatrix;out_els];
load('./data/name_electrodesOnsurface2_0.mat'); % grid 5
elecmatrix=[elecmatrix;out_els];
load('./data/name_electrodesOnsurface3_0.mat'); % grid 6
elecmatrix=[elecmatrix;out_els];
load('./data/name_electrodesOnsurface3_4.mat'); % grid 7
elecmatrix=[elecmatrix;out_els];

% save all projected electrode locaions in a .mat file
save('./data/name_electrodes_surface_loc_all1.mat','elecmatrix');

% make a NIFTI image with all projected electrodes
[output,els,els_ind,outputStruct]=...
    position2reslicedImage2(elecmatrix,'t1_class_gray_electrodes.nii');

for filenummer=1:100
    outputStruct.fname=['electrodes_surface_all' int2str(filenummer) '.img' ];
    if ~exist(outputStruct.fname,'file')>0
        disp(['saving ' outputStruct.fname]);
        % save the data
        spm_write_vol(outputStruct,output);
        break
    end
end

%% 6) generate cortex to render images:
gen_cortex_click('name',0.2,1); 

%% 6) plot electrodes on surface

subject='KB2';
% load cortex
load(['./data/' subject '_cortex.mat']);
% load electrodes on surface
load(['./data/' subject '_electrodes_surface_loc_all.mat']);

ctmr_gauss_plot(cortex,[0 0 0],0)
el_add(elecmatrix,'r',20);
loc_view(270,0)

