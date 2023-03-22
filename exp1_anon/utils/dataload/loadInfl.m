function [inflections]=loadInfl(betaInfo, SIDs, addresp, Dvow, bef, aft)

    % set properties of analysis
    infoStruct = betaInfo; % variable so we can use different speaker betas

    if nargin<3, addresp = 0; end
    if nargin<4, Dvow = []; end

    % ask input from user
    formant = 1;
    % formant to check for no effect electrodes
    formant_noeff = 2;
    rsqthresh = 0; 
    
    % lr_f#     - linear r^2 for F# values
    % sr_f#     - sigmoidal r^2 for F# values
    % infl_f#   - estimated inflection point for F# based on sigmoid fit
    % b_f#      - betas for F# values
    if addresp
        varnames = {'SID', 'el', 'lr_f1', 'sr_f1', 'sl_f1', 'infl_f1', 'b_f1', ...
            'lr_f2', 'sr_f2', 'sl_f2', 'infl_f2', 'b_f2', 'resp'};
    else
        varnames = {'SID', 'el', 'lr_f1', 'sr_f1', 'sl_f1', 'infl_f1', 'b_f1', ...
            'lr_f2', 'sr_f2', 'sl_f2', 'infl_f2', 'b_f2'};       
    end
    inflections =  array2table(zeros(0,length(varnames)), 'VariableNames', varnames);

    for s = 1:length(SIDs)
        SID=SIDs{s};
        
        if isfield(infoStruct, SID) && ~isempty(infoStruct.(SID).els)           
        % difference between sigmoid and linear rsq (for both formants)   
            
            % significant linear vs significant rsq
            lrf1 = squeeze(infoStruct.(SID).rsq(1, :, formant))';
            srf1 = squeeze(infoStruct.(SID).rsq(3, :, formant))';       
            lrf2 = squeeze(infoStruct.(SID).rsq(1, :, formant_noeff))';
            srf2 = squeeze(infoStruct.(SID).rsq(3, :, formant_noeff))';
            
            % slopes
            slf1=infoStruct.(SID).beta_lin(:, formant, 1)>0;
            slf2=infoStruct.(SID).beta_lin(:, formant_noeff, 1)>0;  
            
            % infls
            inf1=infoStruct.(SID).beta_sigm(:, formant, 3);
            inf2=infoStruct.(SID).beta_sigm(:, formant_noeff, 3);
            sids = repmat(str2double(SID(2:end)), length(infoStruct.(SID).els), 1);
            b1=num2cell(squeeze(infoStruct.(SID).y(:, formant, :)), 2);
            b2=num2cell(squeeze(infoStruct.(SID).y(:, formant_noeff, :)), 2);
    
            if addresp
                % neural response for single electrode decoding
                resp = mat2cell(squeeze(Dvow.(SID).resp(infoStruct.(SID).els, :, :)), ...
                    ones(1, length(infoStruct.(SID).els)), bef+aft+1, length(Dvow.vowel));
    
                t2 = table(sids, infoStruct.(SID).els, lrf1, srf1, slf1, inf1, b1, ...
                    lrf2, srf2, slf2, inf2, b2, resp, 'VariableNames', varnames);
            else
                t2 = table(sids, infoStruct.(SID).els, lrf1, srf1, slf1, inf1, b1, ...
                    lrf2, srf2, slf2, inf2, b2, 'VariableNames', varnames);
            end

            inflections = [inflections; t2];
            clear sigrsq* linrsq* 
        end
    end
    
    disp('----------------- Subject Stats -------------------------')
    disp(['Total electrodes: ' num2str(height(inflections))]);
    disp(['Total subj: ' num2str(length(unique(inflections.SID)))]);
    disp(['Min electrodes per subj: ' num2str(min(crosstab(inflections.SID)))]);
    disp(['Max electrodes per subj: ' num2str(max(crosstab(inflections.SID)))]);
    disp(['Median electrodes per subj: ' num2str(median(crosstab(inflections.SID)))]);
end
