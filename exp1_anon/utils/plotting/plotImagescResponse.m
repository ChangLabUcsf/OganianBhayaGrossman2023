%% plot imagesc high gamma response peak across formant space
function plotImagescResponse(DDvow, SID, singleSet, type, ax, cols, bins)

    if nargin<4, type=1; end
    if nargin<6, cols = brewermap(5, 'Dark2'); end
    if nargin<7, bins = []; end
        
    if isempty(bins)
        [N, XEDGES, YEDGES, BINX, BINY] = getF12Bin(DDvow);          
    else
         % clean up formant space
        frm = DDvow.formantVals(1:2, :);
        [~, torem]=rmoutliers(frm(1, :));
        frm(:, torem)=NaN;
        [~, torem]=rmoutliers(frm(2, :));
        frm(:, torem)=NaN;
    
        [N,XEDGES,YEDGES,BINX,BINY]=histcounts2(frm(2, :), frm(1, :), ...
            bins);
    end    
    
    tps = size(DDvow.(SID).resp(:, :, :), 2);
    [~, maxidx]=max(reshape(DDvow.(SID).resp, [], tps), [], 2, 'omitnan');
    peakWind = median(maxidx)-3:median(maxidx)+3;
    
    hgmean=squeeze(mean(DDvow.(SID).resp(:, peakWind, :), 2, 'omitnan'));
    if type==2
        [~, peakIncreaseAdd, ~, ~]=getPeakIncrease(DDvow, 5, SID, DDvow.formantVals, 1, singleSet); 
    end
    
    if nargin < 5
        figure('Position', [100 800 300*length(singleSet) 250]);
        ax = nan(length(singleSet), 1);
    end
    
    ctr=0;    
    for el = singleSet
        ctr=ctr+1;
        img=nan(size(N, 1),size(N, 2));
        for y=1:length(YEDGES)-1
            for x=1:length(XEDGES)-1
                if sum(BINX==x & BINY==y, 'omitnan')>10                
                    switch type
                        case 1
                            % mean hg peak as pixel value
                            img(x, y)=mean(hgmean(el, BINX==x & BINY==y), 'omitnan');
                        case 2
                            % trough to peak as pixel value
                            img(x, y)=mean(peakIncreaseAdd(el, BINX==x & ...
                                BINY==y), 'omitnan');
                    end      
                end
            end
        end
        if nargin < 5
            ax(ctr)=subplot(1, length(singleSet), ctr);
        end
        cla reset;
        imagesc(smoothdata(img')); hold on; 
        
        colormap( [1 1 1; parula(256)]);
               
        ylabel('F1'); xlabel('F2');
        title(num2str(el));  
        
        medoids = plotVariance(DDvow.formantVals, DDvow.vowel, cols, 0, 0); 
        ys=arrayfun(@(x) find(abs(x-YEDGES)==min(abs(x-YEDGES))), medoids(:, 1));
        xs=arrayfun(@(x) find(abs(x-XEDGES)==min(abs(x-XEDGES))), medoids(:, 2));
        scatter(xs, ys, 45, cols, 'filled', 'MarkerEdgeColor', 'k')
        set(gca,'XDir','reverse');
    end
    linkprop(ax, 'clim');
    colormap( [1 1 1; flipud(pink(256))]);   

    % figure formatting
    cbh = colorbar;
    cs = caxis;
    set(cbh,'YTick',[cs(1)+0.05 0 cs(2)-0.05], ...
        'TickLabels', {'min', '0','max'});
    clear cs cbh

    yticks([]); xticks([])
    set(gca, 'FontSize', 15)
    set(gca,'Visible','off');
    
    axes('Position', [ax.Position(1) + 0.01,  ...
        ax.Position(2) .055 .075], 'LineWidth', 2, 'Color', 'none');
    yticks([0 1]); xticks([0 1]);
    xticklabels({'+', '-'});
    yticklabels({'+', '-'});
    ylabel('F1'); xlabel('F2');
    set(gca, 'FontSize', 13)
end