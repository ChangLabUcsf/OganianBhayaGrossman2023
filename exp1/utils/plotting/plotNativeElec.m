%% plot electrodes on NATIVE brains
% desel is like desmat where electrodes are specified under condition
function [native_plot] = plotNativeElec(SIDs, desel, datapath)
    [imgall] = loadElecLoc(datapath);
    
    % last four columns correspond to coordinates and whether electrode is
    % selected (highlighted)
    varnames = {'SID', 'el', 'cond', 'hemi', 'x', 'y', 'z', 'sel', 'anatomy'};
    native_plot = array2table(zeros(0,9), 'VariableNames', varnames);
    for s=1:length(SIDs)
        SID=SIDs{s};  
        
        % check if montage information exists and is on correct hemisphere
        if isfield(imgall, SID)
            img_native=imgall.(SID).img_native;
            
            x_add=20;
            % z dimension for specific subjects    
            if strcmp(SID, 'EC221'), x_add=-90; end
            if strcmp(SID, 'EC100'), x_add=5; end
            if strcmp(SID, 'EC214'), x_add=5; end
            try 
                % plot electrodes corresponding to each condition
                for i=desel.conds
                    
                    % find electrodes corresponding to current condition
                    condel = find(desel.(SID).condition == i);  
                    if ~isempty(condel) && isfield(img_native, 'elecmatrix')
                        
                        elidx = desel.(SID).elid(condel);
                        % make sure not indexing out of mni elecs
                        elidx = elidx(elidx<size(img_native.elecmatrix, 1));
                        n = length(elidx);
                        
                        % plotting in three dimensional space
                        x = img_native.elecmatrix(elidx, 1)+x_add-10;
                        y = img_native.elecmatrix(elidx, 2);
                        z = img_native.elecmatrix(elidx, 3);
                        selid = zeros(n, 1);
                        if isfield(desel.(SID), 'selid')
                            selid = ismember(elidx, desel.(SID).selid);
                        end
                        %x = cat(1, x, [b, c]);                        
                        
                        % make sure dimensions match
                        elidx = reshape(elidx, [n, 1]);
                        selid = reshape(selid, [n, 1]);
                        anatomy = img_native.anatomy(elidx, 4);
                        t2 = table(repmat({SID}, n, 1), elidx, repmat(i, n, 1), ...
                            repmat({imgall.(SID).hemi}, n, 1), x, y, z, selid, ...
                            anatomy, 'VariableNames', varnames);
                        native_plot = [native_plot; t2];   
                        clear x y z n
                    end
                end 
                
                idx = strcmp(SID, native_plot.SID);
                % if there are >0 electrodes to plot for this subject
                if ~isempty(idx)
                    % plot subject specific glass brain
                    fig = figure();
                    axh = axes('Parent', fig);
                    hold(axh, 'all');   
                    ctmr_gauss_plot(img_native.cortex,[0 0 0], 0, imgall.(SID).hemi)

                    condidx = arrayfun(@(x) find(desel.conds==x), native_plot.cond(idx));               
                    scatter3(native_plot.x(idx), native_plot.y(idx), native_plot.z(idx), ...
                        desel.sz(condidx), desel.cols(condidx, :), 'o', 'filled', ...
                        'MarkerFaceAlpha', 0.8); % , 'MarkerEdgeColor', [0 0 0]
                    legend(desel.labels);
                    hold on;

                    grid(axh, 'on');  
                    legend('off');

                    sel = native_plot.sel(idx)>0;
                    scatter3(native_plot.x(sel), native_plot.y(sel), native_plot.z(sel),...
                        75, 'r', 'o', 'MarkerEdgeColor', 'r', 'LineWidth', 1);    
                end
            catch
                warning(['Missing electrode information for ' SID])
            end
        end
    end
               
    disp('------------------- Stats ---------------')
    disp(['Number of subjects included: ' strjoin(unique(native_plot.SID))]);
    disp(['Total number of electrodes: ' num2str(height(native_plot))]);
    disp(['Number of electrodes over cond 1: ' num2str(sum(native_plot.cond>1))]);   
    
end