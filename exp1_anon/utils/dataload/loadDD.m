%% load DD info
function Dvow = loadDD(corpus, bef, aft, SIDs, datapath)

    addpath(genpath(datapath));
    alloudness = getfield(load(['stim_info/out_sentence_details_' corpus ...
        '_all_loudness.mat']), 'sentdet'); 
    Dvow = load(['stim_info/Dvow_' corpus]);
    
    % populates formant values
    Dvow.formantVals = nan(4, length(Dvow.vowel));    
    for sent=1:length(alloudness)
        sentName=alloudness(sent).name;
        sentIdx=find(ismember({alloudness.name},sentName));
        vowId=alloudness(sentIdx).vowelId;   
        Dvow.formantVals(:, vowId) = alloudness(sentIdx).frmMedVal;  
    end
    
    % populates mean pitch values
    fid = fopen(['stim_info/info_meanF0_' corpus '.txt']);
    a = textscan(fid, '%s%f', 'Delimiter', ':');
    fclose(fid); meanf0.name = a{1}; meanf0.f0 = a{2};
    for sent=1:length(meanf0.name)
        [mem, sentId]=ismember(meanf0.name{sent}(1:end-1), {alloudness(:).name});
        if mem
            Dvow.meanf0(1, Dvow.sent==sentId) = ... dddd
                repmat(meanf0.f0(sent), 1, sum(Dvow.sent==sentId, 'omitnan'));
        end
    end    
    
    % add all subject specific resp matrices to Dvow  
    Dvow.repVows = zeros(length(Dvow.vowel), 1);
    Dvow=addtoDD(Dvow, corpus, bef, aft, SIDs, alloudness, 0, datapath);
    
    % adds normalized formant values to Dvow
    [Dvow.normformantVals, ~] = normFormantVals(Dvow, corpus, 1, 'formant0', 2); 
    
    % adds vowelType to Dvow
    Dvow.vowelType = zeros(1, length(Dvow.vowel));
    for n = unique(Dvow.vowel)
        vowIdx = find(ismember(unique(Dvow.vowel), n));
        Dvow.vowelType = vowIdx*ismember(Dvow.vowel, n)+Dvow.vowelType;
    end

end