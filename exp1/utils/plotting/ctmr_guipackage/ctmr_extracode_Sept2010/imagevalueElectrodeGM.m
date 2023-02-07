function output=imagevalueElectrodeGM(circleradius,els,tName,gmName)
% function output=imagevalueElectrode(circleradius)
% function calculated average value of map (e.g.Tmap) withing circle of circleradius for
% native space coordinates els
%
% input:
% circleradius=5;
%
% output:
% output.mean - mean value per electrode
% output.max - max value per electrode

% gray matter threshold
gmThreshold=0.1;

%select T map
%data.tName=spm_select(1,'image','select map for values (NONresliced t-map/percsignalchange)');
data.tStruct=spm_vol(tName);
[data.t,data.txyz]=spm_read_vols(data.tStruct);% from structure to data matrix
data.txyz=data.txyz';

data.gmStruct=spm_vol(gmName);
[data.gm]=spm_read_vols(data.gmStruct);

% clear T-values not in gray matter:
data.t(data.gm<gmThreshold)=NaN;

output.mean=zeros(size(els(:,1)))+7;
output.max=zeros(size(els(:,1)))+7;
output.min=zeros(size(els(:,1)))+7;

%for each electrode draw circle and get average T-value
for elec=1:length(els)
    temp.electrode=zeros(size(data.t));
    temp.dist=sqrt((data.txyz(:,1)-els(elec,1)).^2+...
        (data.txyz(:,2)-els(elec,2)).^2+...
        (data.txyz(:,3)-els(elec,3)).^2);
    temp.electrode(temp.dist<circleradius)=1;
    temp.T=data.t(temp.electrode>0);
    temp.T=temp.T(~isnan(temp.T));
    if ~isempty(temp.T)
        output.mean(elec)=mean(temp.T);
        output.max(elec)=max(temp.T);
        output.min(elec)=min(temp.T);
    else
        output.mean(elec)=NaN;
        output.max(elec)=NaN;
        output.min(elec)=NaN;
    end
end

