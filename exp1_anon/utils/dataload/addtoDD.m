%% add subject information to DD struct
function Dvow = addtoDD(Dvow, corpus, bef, aft, SIDs, alloudness, repeats, datapath)
    if nargin<8, repeats = 0; end 
    
    % check if SIDs have already been added    
    toAdd = cellfun(@(x) ~isfield(Dvow, x), SIDs);
    
    if ~any(toAdd)
        warning('All SIDs already added.');
        return
    end
    
    if nargin<7 || isempty(alloudness), alloudness = ...
            getfield(load(['out_sentence_details_' corpus '_all_loudness.mat']),'sentdet'); end
    
    for subj=find(toAdd)
        disp(['loading ...... subj ' num2str(subj) '/' num2str(length(SIDs))])
        
        % additional check if SID already exists in current DVow struct
        if ~isfield(Dvow, SIDs{subj})            
            %%  load subject data            
            
            tic
            if ~repeats                              
                % align and fill DDvow
                filename = fullfile(datapath, 'pt_data', SIDs{subj}, ...
                    'dimex', [SIDs{subj} '_Dvow.mat']);
                if isfile(filename)
                    load(filename, 'resp');
                else % only works with access to full sentence responses 
                    out = addNeuralResp(SIDs{subj}, corpus, datapath);
                    
                    els = size(out(1).resp, 1);
                    tps = bef + aft + 1;

                    Dvow.(SIDs{subj}).resp = nan(els, tps, length(Dvow.vowel));
                    Dvow.(SIDs{subj}).predResp = nan(els, tps, length(Dvow.vowel));
        
                    for sent=1:length(out)
                        sentName=out(sent).name;
                        sentIdx=find(ismember({alloudness.name},sentName));
                        
                        vowIdx=alloudness(sentIdx).vowelTimes(1, :);
                        vowId=alloudness(sentIdx).vowelId;  
    
                        % previous check that was a bug -if vowIdx(vow)-bef>0 && ...
                        % vowIdx(vow)+aft < size(out(sent).resp, 2)+1 
                        outresp = arrayfun(@(x) replaceMat(out(sent).resp(:, ...
                            vowIdx(x)-bef:vowIdx(x)+aft), 0, NaN), 1:length(vowId), ...
                            'UniformOutput',false);
                        resp(:, :, vowId)= cat(3, outresp{:});                                                                              
                    end
                    save(filename, 'resp')
                end

                % pad array with nans if not all vowels were presented
                if length(Dvow.vowel)-size(resp, 3)>0
                    Dvow.(SIDs{subj}).resp = padarray(resp, ...
                        [0, 0, length(Dvow.vowel)-size(resp, 3)], NaN, 'post');
                else
                    Dvow.(SIDs{subj}).resp = resp;
                end
                Dvow.(SIDs{subj}).resp(Dvow.(SIDs{subj}).resp==0)=NaN;
            else

                % set up for aligned responses
                out = addNeuralResp(SIDs{subj}, corpus, datapath);

                if ~isempty(out)
                    els = size(out(1).resp, 1);
                    tps = bef + aft + 1;
    
                    repSent  = find(cellfun(@(x) length(x)>9, {out.Trials}));
                    repNum = cellfun(@(x) length(x), {out(repSent).Trials});
                    repNum = min(repNum);
                    Dvow.(SIDs{subj}).resp = nan(els, tps, length(Dvow.vowel), repNum);
                    Dvow.(SIDs{subj}).repVows = zeros(1, length(Dvow.vowel));
    
                    % align and fill DDvow
                    for sent=repSent
                        sentName=out(sent).name;
                        sentIdx=find(ismember({alloudness.name},sentName));
                        
                        vowIdx=alloudness(sentIdx).vowelTimes(1, :);
                        vowId=alloudness(sentIdx).vowelId;   
                        
                        outresp = arrayfun(@(x) replaceMat(out(sent).resp(:, ...
                            vowIdx(x)-bef:vowIdx(x)+aft, 1:repNum), 0, NaN), ...
                            1:length(vowId),'UniformOutput',false);
                        Dvow.(SIDs{subj}).resp(:, :, vowId, :)= permute(...
                            cat(4, outresp{:}), [1, 2, 4, 3]);
                        clear outresp                                                                                                            
                    end
                end
            end
            toc
        end     
    end
end