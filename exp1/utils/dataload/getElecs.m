function [allidx, fvals] = getElecs(data, SIDs, corpus, type, ...
    typeinfo, datapath)
% saved out to 'out_elecs_....'
% function [allidx] = get_elecs()
% selects electrode based on the typeinfo structure
% Inputs: 
            % data - out structure, contains a 'resp' field under SID field with neural
            % data. resp is of form electrodes x timepoints x trials
            % SIDs - SIDs to run selection on
            % tps - time window for ftest (will be averaged over) or
            % 'bysubj', 'bychan'
            % corpus - 'dimex' or 'timit'
            % type - 'ftest', 'anatomy', 'strfrsq'
            % typeinfo    
            %     'ftest' - contains sigField region with labels for
            %     anovan, and fieldname for response (resp)
            %     'anatomy' - contains region field with anatomical roi
            %     'strfrsq' - contains modelname and rsqthresh values for
            %     strf and rsq threshold value
% Outputs:
            % allidx - struct with SID fields containing electrode indices
            % of selected electrodes
            % fvals - if type is ftest, returns all electrode F-values that
            % were thresholded to find subset of selected electrodes
% Ilina Bhaya-Grossman, August 2020
    
    % how should electrodes be selected? 
    switch type                    
        case 'strfrsq'
            % uses model strf > rsqthresh
            if ~isfield(typeinfo, 'zscoreflag')
                typeinfo.zscoreflag = 1;
                warning('Default zscore flag will be used.');
            end
            
            if ~isfield(typeinfo, 'modelname')
                error('Please add a strf modelname to the input struct.');
            end
            
            if ~isfield(typeinfo, 'rsqthresh')
                typeinfo.rsqthresh=0.1;
                warning('Default rsq threshold will be used.');
            end
        case 'anatomy'
            if ~isfield(typeinfo, 'region')
                typeinfo.region='superiortemporal';
                warning('Default region will be used.');
            end            
            typeinfo.imgdata = loadElecLoc(datapath);
    end
    
    if ~isfield(typeinfo, 'resp')
        typeinfo.resp = 'resp';
        warning('Default response field will be used.');
    end
    allidx = [];
    fvals = [];
       
    ctr=1;
    % turn on to see Fstat plots
    % loop over SIDs to aggregate all electrodes that meet type requirements
    for subj=1:length(SIDs)
       SID=SIDs{subj};   
       disp(['loading electrodes ...... subj ' num2str(subj) '/' num2str(length(SIDs))])
        
        switch type
            case 'strfrsq'                
                strf=loadMultModelStrf(SID, typeinfo.modelname, corpus, ...
                    datapath, typeinfo.zscoreflag, '');
                if any(cellfun(@(x) isempty(x), strf))
                    warning(['No strf model run for ' SID]);
                    elidx = [];
                else                   
                    switch length(typeinfo.modelname)                
                        case 1 % all el where model 1 > thresh
                            elidx=find(strf{1}.meanTestR.^2>typeinfo.rsqthresh);
                        case 2 % all el where model 2 > model 1 and model 2 > thresh
                            elidx = find(strf{2}.meanTestR.^2 > typeinfo.rsqthresh ...
                                & strf{2}.meanTestR>strf{1}.meanTestR);  
                        case 4 % all el where model 2 > model 1, model 4 > model 3 and model 2 > thresh       
                            elidx = find(strf{2}.meanTestR.^2 > typeinfo.rsqthresh ...
                                & strf{2}.meanTestR>strf{1}.meanTestR ...
                                & strf{4}.meanTestR>strf{3}.meanTestR); 
                    end
                    % second argument returned is the first model rsqs
                    fvals.(SID) = strf{1}.meanTestR.^2;
                end
            case 'anatomy'
                elidx=find(strcmp(typeinfo.imgdata.(SID).imgNative.anatomy(:, 4), ...
                    typeinfo.region));
                % in case there are stereo electrodes we are not including
                elidx=elidx(elidx<size(data.(SID).resp, 1)); 
        end    
        
        if ~isempty(elidx)                       
            % add electrode and SID information to allidx array
            allidx.(SID) = elidx;
        else
            warning(['No electrodes from ' SID ' will be included.'])
        end                
        ctr=ctr+length(elidx);       
        
        clear pval elidx SID
    end
    
end
