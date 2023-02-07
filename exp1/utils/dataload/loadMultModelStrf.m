%% get corpusStrfs
function corpusStrf=loadMultModelStrf(SID, modelnames, corpus, datapath, zscoreflags, pitchSel)

    if nargin<5, zscoreflags=0; end

    % pitch selection criteria for sentences included in training + test
    % set
    if nargin<7, pitchSel=''; end % 'lte170', 'gt170', ''
    
    corpusStrf = cell(1, length(modelnames));
    for i=1:length(modelnames)
        if length(zscoreflags)>1 % can use different zscore flags for the multiple models
            zscoreflag = zscoreflags(i);
        else
            zscoreflag = zscoreflags;
        end
        scaleflag=0; % version is v5

        modelname = modelnames{i};
        
        pitchStr = pitchSel;
        if ~strcmp(pitchSel, '')
            pitchStr = ['_' pitchSel]; 
        end
        
        filename = [SID '_strf_' modelname '_' ...
            'zX' num2str(zscoreflag) '_zY0_scX1_scY0_hp0_SentOns1_sentScale' ...
            num2str(scaleflag) '_El1to*_600ms_edge1_boot1' pitchStr '.mat'];
        fulldir = dir(fullfile(datapath, SID, [corpus '/strf/'], filename));
        
        try
            corpusStrf{i} = load(fullfile(fulldir.folder, fulldir.name));
        catch
            warning(['No ' corpus ' ' modelname ' strf data exists for ' SID])
        end
    end
end


