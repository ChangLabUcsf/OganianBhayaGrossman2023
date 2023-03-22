%% normalize formant medians by sentence or f0 (sentence is a subset of f0)
function [formantValues, Y]=normFormantVals(Dvow, corpus, formant, normBy, bins)
    % normalize formant values by sentence
    Y=[];
    if nargin<5, bins = 1; end
    if nargin<4, normBy = 'sentence'; end
    formantValues = nan(size(Dvow.formantVals));
    
    if strcmp(normBy, 'sentence')
        sentdet = load(['out_sentence_details_' corpus '_all_loudness.mat']);  
        sentdet = sentdet.sentdet;
        for s = 1:length(sentdet)
            formNan = sentdet(s).formants;
            formNan(formNan==0) = NaN;
            vowelTimes = sentdet(s).vowelTimes;

            % norm over vowel time points
            allVowelTimes = [];
            for v = 1:length(vowelTimes)
                tmp = vowelTimes(1, v):vowelTimes(2, v);
                allVowelTimes = [allVowelTimes tmp];
            end

            % norm over entire time course
            formNan(:, ~ismember(1:size(formNan,2), allVowelTimes)) = NaN;
            normFrm = normalize(formNan, 2);        

            vowelId = sentdet(s).vowelId;
            for v = 1:length(vowelTimes)
                tmp = normFrm(formant, vowelTimes(1, v):vowelTimes(2, v));
                formantValues(vowelId(v)) = median(tmp);
            end
        end
    else
        % trying this without the entire vowel time course
        Y = get_densbin_arr(Dvow.meanf0,Dvow.meanf0,bins); 
        for s = 1:bins 
            % norm over entire time course
            normFrm = normalize(Dvow.formantVals(:, Y==s), 2);  
            % normalizing over all formants (fixed from previous when one formant was normed)
            formantValues(:, Y==s) = normFrm;
        end
    end
end
