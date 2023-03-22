function out_els=p_zoom_array(els, gs, lcl_num)
% function out_els=p_zoom(els, gs, lcl_num);
% this function projects electrodes to hemispheric surface hull "gs" (verts x 3) that is
% closest for each electrode in "els" (N x 3), in direction of local norm
% lcl_num indicates the number of point for calculation of norm, 0 if global
%   Created by:
%   D. Hermes & K.J. Miller 
%   Dept of Neurology and Neurosurgery, University Medical Center Utrecht
%
%   Version 1.1.0, released 26-11-2009
%   modified by kjm 4/2011

%% get point at center of hemispheric hull
    cent_pt=mean(gs,1); %center point
    gs_dist_to_cent = sum((gs-ones(size(gs,1),1)*cent_pt).^2,2).^.5;

%% get global norm vector
    if lcl_num==0 %global estimate of principal direction most orthogonal to array
        [v,d]=eig(cov(els)); %all vecs
        [y,ind]=min(abs(diag(d))); nm=v(:,ind); %vec we want
        mg = mean(els,1); %mean pt in grid
        if dot(nm,(mg-cent_pt))<0, nm = -nm; end % make sure that vector points away from the center of the hull
    end

%%
out_ind=zeros(size(els,1),1);

for k=1:size(els,1)
    
    %% sub array? get local norm vector
    if lcl_num>0, % get principal direction most orthogonal to sub-array
        [y,ind]=sort(dist(els,els(k,:).'),'ascend');%select closest for sub-array
        els2use=ind(2:(lcl_num+1)); %use lcl number of points (note that we exclude the important one)
        %
        [v,d]=eig(cov(els(els2use,:))); %all vecs,
        [y,ind]=min(abs(diag(d))); nm=v(:,ind); %vec we want
        if dot(nm,(els(k,:)-cent_pt))<0, nm = -nm; end % make sure that vector points away from the center of the hull   
    end
    
    %%
    % get unit vector in direction to every point on hull
    npls=[gs(:,1)-els(k,1) gs(:,2)-els(k,2) gs(:,3)-els(k,3)]; %x,y,z vector from each point to our point of interest
    tmp=(sum(npls.^2,2).^.5);
    npdot=npls./(tmp*[1 1 1]);
    % project unit vectors onto norm vector
    tmp=dot(npdot,ones(size(gs,1),1)*(nm.'),2);
    
    
    if max(tmp<.65)  % fix for case that our electrode is outside of the hull to begin with (which shouldn't typically happen) 
         % limit distance to electrode
        gs_dist_to_el = sum((gs-ones(size(gs,1),1)*els(k,:)).^2,2).^.5;
        tmp=abs(tmp).*(1./gs_dist_to_el);
    end
        
        % get projected to point
        [y,out_ind(k)]=max(tmp);
    
end

out_els=gs(out_ind,:);


% plot on surface to check
figure
plot3(els(:,1),els(:,2),els(:,3),'b.','MarkerSize',20);
hold on;
plot3(out_els(:,1),out_els(:,2),out_els(:,3),'r.','MarkerSize',20);
plot3(gs(:,1),gs(:,2),gs(:,3),'k.','MarkerSize',1);
