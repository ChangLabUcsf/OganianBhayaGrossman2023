function [xyztCortex,t_surf] = voxplot_func_gm(Tthreshold,Dthreshold,varargin)
%
% use as follows:
% [xyztCortex,t_surf] = voxplot_func_gm(Tthreshold,Dthreshold,[]);
% or
% [xyztCortex,t_surf] = voxplot_func_gm(Tthreshold,Dthreshold,GMthreshold);
%
% input
% Tthreshold - threshold for t statistic
% Dthreshold - depth under surface
%
% optional: GM threshold - if specified, only plots voxels from t-map in gray matter
% >= GMthreshold
%
% output used in ctmr_vox_plot to make a rendering
%
%      Copyright (C) 2009  D. Hermes, Dept of Neurology and Neurosurgery, University Medical Center Utrecht
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


%% load files

% load surface
[sName]=spm_select(1,'image','select image with surface');
% load tmap
[tName]=spm_select(1,'image','select image with tmap');

if ~ isempty(varargin{1})
    GMthreshold=varargin{1};
    % ask for gray matter
    [gmName]=spm_select(1,'image','select image with rc1');
end

% load cortex
[cname]=spm_select(1,'mat','select cortex.mat');

load(cname); % cortex
s_info=spm_vol(sName); [s]=spm_read_vols(s_info); % surface
t_info=spm_vol(tName); [t]=spm_read_vols(t_info); % t statistics

if ~ isempty(varargin{1})
    gm_info=spm_vol(gmName); [gm]=spm_read_vols(gm_info);
    % clear T-values not in gray matter:
    t(gm<GMthreshold)=NaN;
end



%%
% indices to native surface
[xs,ys,zs]=ind2sub(size(s),find(s>0));
xyzs=[xs ys zs];
clear xs ys zs % housekeeping
xyzs=(xyzs*s_info.mat(1:3,1:3)')+repmat(s_info.mat(1:3,4),1,length(xyzs))';

% indices to native tmap
[xt,yt,zt]=ind2sub(size(t),find(t>=Tthreshold));
xyzt=[xt yt zt];
clear xt yt zt % housekeeping
xyzt=(xyzt*t_info.mat(1:3,1:3)')+repmat(t_info.mat(1:3,4),1,length(xyzt))';

% select tmap coordinates that are within Dthreshold mm of surface
tsel=zeros(length(xyzt),1);
tnewind=zeros(length(xyzt),1);
for k=1:length(xyzt)
    distvect=(sum((xyzs-repmat(xyzt(k,:),length(xyzs),1)).^2,2)).^0.5;
    if min(distvect)<=Dthreshold
        tsel(k)=1;
        tnewind(k)=find(distvect==min(distvect),1);
    end
end
clear distvect;
tnewind=tnewind(tnewind>0);

% set rest to []
xyzt(tsel==0,:)=[];
xyztSurf=xyzs(tnewind,:);

% get t values
t_surf=t(t>=Tthreshold);
t_surf=t_surf(tsel>0);

% get coordinates on cortex (instead of surface)
xyztCortex_ind=zeros(length(xyzt(:,1)),1);
for k=1:length(xyzt(:,1))    
    distvect=(sum((cortex.vert-repmat(xyzt(k,:),length(cortex.vert),1)).^2,2)).^0.5;
    xyztCortex_ind(k)=find(distvect==min(distvect),1);
end
xyztCortex=cortex.vert(xyztCortex_ind,:);

