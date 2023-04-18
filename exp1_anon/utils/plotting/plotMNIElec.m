%% plot electrodes on MNI brain
% desel is like desmat where electrodes are specified under condition
function [mni_plot] = plotMNIElec(SIDs, desel, hemi, datapath)
    [imgall] = loadElecLoc(datapath);

    fig = figure();
    axh = axes('Parent', fig);
    hold(axh, 'all');
        
    SIDs = SIDs(cellfun(@(x) strcmp(imgall.(x).hemi, hemi), SIDs));  
    if isempty(SIDs)
        warning('No SIDs included with correct hemisphere type.');
        mni_plot = [];
        return
    end 

    % plot glass brain    
    ctmr_gauss_plot(imgall.(SIDs{1}).img_mni.cortex,[0 0 0], 0, hemi)

    % last four columns correspond to coordinates and whether electrode is
    % selected (highlighted)
    varnames = {'SID', 'el', 'cond', 'hemi', 'x', 'y', 'z', 'sel', 'anatomy'};
    mni_plot = array2table(zeros(0,9), 'VariableNames', varnames);
    for s=1:length(SIDs)
        SID=SIDs{s};
        
        % check if montage information exists
        if isfield(imgall, SID)
            img_mni=imgall.(SID).img_mni;
            
            x_add=70;
            % z dimension for specific subjects
            if ismember(SID, {'S6', 'S5'}), x_add=40; end
            if ismember(SID, {'S1'}), x_add=-50; end
            if ismember(SID, {'S9'}), x_add=-200; end
            if strcmp(hemi,'lh'), x_add=x_add-90; end
            
            try 
                % plot electrodes corresponding to each condition
                for i=desel.conds
                    
                    % find electrodes corresponding to current condition
                    condel = find(desel.(SID).condition == i);  
                    if ~isempty(condel)
                        
                        % electrodes in condition, make sure not indexing out of mni elecs
                        elidx = desel.(SID).elid(condel);                        
                        elidx = elidx(elidx<size(img_mni.elecmatrix, 1));
                        n = length(elidx);
                        
                        % plotting in three dimensional space
                        x = img_mni.elecmatrix(elidx, 1)+x_add-10;
                        y = img_mni.elecmatrix(elidx, 2);
                        z = img_mni.elecmatrix(elidx, 3);
                        selid = zeros(n, 1);
                        if isfield(desel.(SID), 'selid')
                            selid = ismember(elidx, desel.(SID).selid);
                        end                      
                        
                        % make sure dimensions match
                        elidx = reshape(elidx, [n, 1]);
                        selid = reshape(selid, [n, 1]);
                        anatomy = img_mni.anatomy(elidx, 4);
                        t2 = table(repmat({SID}, n, 1), elidx, repmat(i, n, 1), ...
                            repmat({hemi}, n, 1), x, y, z, selid, anatomy, 'VariableNames', varnames);
                        mni_plot = [mni_plot; t2];   
                        clear x y z n
                    end
                end 
            catch
                warning(['Missing electrode information for ' SID])
            end
        end
    end
    
    % convert the conditions to be indices for the size field
    % unique sorts
    condidx = arrayfun(@(x) find(unique(desel.conds)==x), mni_plot.cond);    
    scatter3(mni_plot.x, mni_plot.y, mni_plot.z, desel.sz(condidx), ...
        desel.cols(condidx, :), 'o', 'filled', ...
        'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0 0 0]);
  
    % for selected electrodes
    sel = mni_plot.sel>0;
    scatter3(mni_plot.x(sel), mni_plot.y(sel), mni_plot.z(sel),...
        75, 'r', 'o', 'MarkerEdgeColor', 'r', 'LineWidth', 1);                    

    if ~isempty(mni_plot)
        disp('------------------- Stats ---------------')
        disp(['Number of subjects included: ' strjoin(unique(mni_plot.SID))]);
        disp(['Total number of electrodes: ' num2str(height(mni_plot))]);
        disp(['Number of electrodes over cond 1: ' num2str(sum(mni_plot.cond>1))]);   
    end 
end