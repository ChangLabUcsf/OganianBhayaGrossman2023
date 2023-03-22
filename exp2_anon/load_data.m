%% load synthetic stimulus details
taskstims = getfield(load('stimdet_vowels_v5.mat'), 'stimdet');
%% load dimex stimulus details
outD.sent = getfield(load(fullfile(stimFolder, 'out_sentence_details_dimex_all_loudness.mat')),'sentdet');

%% task stimuli organization
%  remove high frequencies from aud
for i = 1:length(taskstims)
    taskstims(i).aud(50:end,:) = 0;
end

% create single timepoint aud
for i = 1:length(taskstims)
    taskstims(i).aud1tp = zeros(size(taskstims(i).aud));
    taskstims(i).aud1tp(:,51) = mean(taskstims(i).aud(:,51:end-50),2);    
end

%% aud center frequencies
centerfs = load(fullfile(stimFolder, 'mel_centerF.mat'));
centerfs.meanfs = centerfs.binfrqs(2:end-1) + diff(centerfs.binfrqs(2:end))/2;
%% load dimex vowel info
Dvow = load(fullfile(stimFolder, 'Dvow_dimex.mat'));
Dvow.frm = [outD.sent.frmMedVal];
%% get median vowel formants
for i = 1:4
    for cv = 1:5
    Dvow.frmMedian(i,cv) = median(Dvow.frm(i,strcmpi(Dvow.vowel, vowNames{cv})));
    end
end

%% ---- load segmented data 
for si = 1:length(SID)    
    fprintf(2, '%d: %s ...', si, SID{si});
    DDfilename = [SID{si} '_' taskname '_DD.mat'];
    ddfile = fullfile(outpath, SID{si},DDfilename);  
    alldat(si) = load(ddfile);  
end
%% add stimId to alldat
mistr = [];
for cs = 1:length(alldat)
    for ct = 1:length(alldat(cs).mean.f1)
        ctr = find([taskstims.f1] == alldat(cs).mean.f1(ct) & [taskstims.f2] == alldat(cs).mean.f2(ct));
        if ~isempty(ctr)
            alldat(cs).mean.stimId(ct) = ctr;
        else
            alldat(cs).mean.stimId(ct)=nan;
           mistr(end+1,:) =  [cs, alldat(cs).mean.f1(ct), alldat(cs).mean.f2(ct)];
        end        
    end
end


%% mark points inside/outside contour in alldat
for si=1:length(alldat)
    stims = [alldat(si).mean.f2; alldat(si).mean.f1];
    pgon = polyshape(dimCont.M(1,2:end)', dimCont.M(2, 2:end)');
    alldat(si).mean.inVowSpace = isinterior(pgon, stims(1,:)', stims(2,:)');
end

