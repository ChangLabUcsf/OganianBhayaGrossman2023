
% bins by density (every bin should have the same freq count)
function [binVal, ranges, nanIdx] = get_densbin_arr(valuesToBin, valuesToBinBy, numBins)    
    
    nanIdx = isnan(valuesToBin);
    allVal = sort(valuesToBinBy); 
    allVal(isnan(allVal)) = [];
    
    bin_sz = floor(length(allVal)/numBins);
    binVal = nan(1, size(valuesToBin, 2));

    ranges = {};
    %create formant binned binary features
    for b = 1:numBins             
        % create bins (clean up)
        startEl=-Inf; endEl=Inf; 
        if b~=1, startEl=allVal(bin_sz*(b-1)); end
        if b~=numBins, endEl=allVal(bin_sz*b); end
        ranges(b) = {num2str(round([startEl endEl]))};
%         ranges(b) = {num2str(round([startEl endEl], 2,'significant'))};
        binVal(1, find(inrange(valuesToBin, [startEl endEl]))) = b;
    end 
end