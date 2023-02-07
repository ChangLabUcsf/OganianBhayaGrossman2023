%% plot imagesc average high gamma response peak across formant space
function avgVowelVal = plotAvgImagescResponse(Dvow, allidx, cols, bins, weights)

    if nargin<3, cols = brewermap(5, 'Dark2'); end
    if nargin<4, bins = []; end
    if nargin<5, weights = ones(height(allidx), 1); end

    avgVowelVal = nan(5, 1);
        
    if isempty(bins)
        [~, XEDGES, YEDGES, BINX, BINY] = getF12Bin(Dvow);          
    else
         % clean up formant space
        frm = Dvow.formantVals(1:2, :);
        [~, torem]=rmoutliers(frm(1, :));
        frm(:, torem)=NaN;
        [~, torem]=rmoutliers(frm(2, :));
        frm(:, torem)=NaN;
    
        [~,XEDGES,YEDGES,BINX,BINY]=histcounts2(frm(2, :), frm(1, :), ...
            bins);
    end    
    
    img = nan(height(allidx), length(XEDGES)-1, length(YEDGES)-1);
    for e = 1:height(allidx) 
        SID = allidx{e, 'SID'}{1};
        el = allidx{e, 'el'};
        tps = size(Dvow.(SID).resp, 2);
        [~, maxidx]=max(reshape(Dvow.(SID).resp, [], tps), [], 2, 'omitnan');
        peakWind = median(maxidx)-2:median(maxidx)+2;

        hgmean=squeeze(mean(Dvow.(SID).resp(:, peakWind, :), 2, 'omitnan')); 
        for y=1:length(YEDGES)-1
            for x=1:length(XEDGES)-1
                if sum(BINX==x & BINY==y, 'omitnan')>10                
                    % mean hg peak as pixel value
                    img(e, x, y)=mean(hgmean(el, BINX==x & ...
                        BINY==y), 'omitnan')*weights(e);     
                end
            end
        end              
    end
    
    A = smoothdata(squeeze(mean(img, 1)), 'SmoothingFactor',0.5)';
    imagesc(XEDGES, YEDGES, A); hold on;
    colormap( [1 1 1; parula(256)]);
    ylabel('F1'); 
    xlabel('F2'); 

    % quantile is vowel x quantile x formant
    [medoids, quant] = plotVariance(Dvow.formantVals, Dvow.vowel, cols, 1, 0); 

    for v = 1:5
        center = fliplr(medoids(v, 1:2)); 
        axes = mean(abs(center'-flipud(squeeze([quant(v, [1, 6], :)])')), 2);
        roi = drawellipse('Center',center,'SemiAxes',axes','Color',cols(v, :), ...
            'LineWidth',0.1);
        roi.Visible = 'off';

        avgVowelVal(v) = mean(A(createMask(roi)));
    end
    
    scatter(medoids(:, 2), medoids(:, 1), 45, cols, 'filled', 'MarkerEdgeColor', 'k')
    set(gca,'XDir','reverse');
    colormap( [1 1 1; flipud(pink(256))]);
end