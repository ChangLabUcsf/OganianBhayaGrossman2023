%% plot erps over vowels
function plotVowelErp(Dvow, SID, elect, f, cols, bef)
    
    desmat.names = [];
    desmat.names = unique(Dvow.vowel);
    desmat.condition = Dvow.vowelType;
    unique_vowels = {'a', 'e', 'i', 'o', 'u'};
      
    vowelNameIdx=ismember(unique(Dvow.vowel), unique_vowels);
    names = unique(Dvow.vowel);
    desmat.names = names(vowelNameIdx);

    % all indices with vowels to keep
    vowelCorpusIdx = ismember(Dvow.vowelType,find(vowelNameIdx));

    % relabel vowels in order of how they appear in vowel names
    desmat.condition = nan(1, sum(vowelCorpusIdx, 'omitnan'));
    for i=1:length(find(vowelNameIdx))
        vowelInds = find(vowelNameIdx);
        desmat.condition(1, Dvow.vowelType(vowelCorpusIdx)==vowelInds(i))=i;
    end
    desmat.vowelNames = Dvow.vowel(vowelCorpusIdx);

    resp = Dvow.(SID).resp(:, :, vowelCorpusIdx);
    ecog_erpPlotSingle(resp, desmat, elect, elect,[], f, bef, [], 1, 100, [], cols);   
end