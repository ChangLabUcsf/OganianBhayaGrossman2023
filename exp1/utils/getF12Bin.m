function [N, XEDGES, YEDGES, BINX, BINY] = getF12Bin(DDvow)

    % all formants
    frm = DDvow.formantVals(1:2, :);

    % clean up formant space
    [~, torem]=rmoutliers(frm(1, :));
    frm(:, torem)=NaN;
    [~, torem]=rmoutliers(frm(2, :));
    frm(:, torem)=NaN;

    % use histcounts2 to bin the space
    [N,XEDGES,YEDGES,BINX,BINY]=histcounts2(frm(2, :), frm(1, :), ...
        'BinMethod', 'fd');
end