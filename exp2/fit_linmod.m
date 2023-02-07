%% ------ calculate Linear models with integration across bands
%% mean resp
for csubj = 1:length(alldat)
    alldat(csubj).mean.meanresp = squeeze(nanmean(alldat(csubj).mean.resp(:,alldat(csubj).befaft(1)+(10:40),:),2));
end

%% fit models on mean responses
allprednames = {'Intercept', 'F1', 'F2', 'IA'};
linmDefine = {'F1&F2', 1:3 ; ...
    'F1&F2&IA', 1:4;...
    'intercept', 1;...
    'F1',[1 2];...
    'F2',[1 3]};

for csubj = 1:length(alldat)
    cdat = alldat(csubj);
    %% hga - y
    y = squeeze(nanmean(cdat.mean.resp(:,cdat.befaft(1)+(10:40),:),2));
    cy = y';
    incltr = isfinite(cy(:,1));
    cy(incltr,:) = zscore(cy(incltr,:),0,1);
    %% all predictors 
    allpred = nan(length(cdat.mean.f1), length(allprednames));
    allpred(:,1) = ones(length(cdat.mean.f1), 1);
    allpred(:,2) = zscore(cdat.mean.f1);
    allpred(:,3) = zscore(cdat.mean.f2);
    allpred(:,4) = allpred(:,2).*allpred(:,3);
%     
    %% all models, predictor selection
    for cm  = 1:size(linmDefine,1)
        cpred = allpred(:,linmDefine{cm,2});
        curmfit.bnames = allprednames(linmDefine{cm,2});
        curmfit.b = nan(size(y,1), size(cpred,2),size(y,2));
        curmfit.rsq = nan(cdat.elcount,1);
        curmfit.RSS = nan(cdat.elcount,1);
        curmfit.y_pred = nan(size(cy));
        curmfit.modelName = linmDefine{cm,1};
        % cross-validated parameters
        [curmfit.y_pred(incltr,:),curmfit.b(:,:,incltr), curmfit.rsq, curmfit.RSS,curmfit.rCorr,curmfit.Ftest, curmfit.redM] =...
            crossvalFit_v3_vect(cpred(incltr,:),cy(incltr,:), 'linear');
        % write it out
        linmfits(csubj,cm) = curmfit;
        clear curmfit;
    end
    
end

%% fit model with predictor for in/out of vowel space
vowspprednames = {'Intercept', 'F1', 'F2', 'IA', 'vowSP', 'vowSpxF1', 'vowSpxF2', 'vowSpxIA'};
linmDefine2 = {'(F1+F2+IA)*vowSp', 1:8;...
    'F1+F2+IA+vowSp', 1:5};
cmnum = size(linmfits,2);
for csubj = 1:length(alldat)
    cdat = alldat(csubj);
    %% hga - y
    y = squeeze(nanmean(cdat.mean.resp(:,cdat.befaft(1)+(10:40),:),2));
    cy = y';
    incltr = isfinite(cy(:,1));
    cy(incltr,:) = zscore(cy(incltr,:),0,1);
    %% all predictors 
    allpred = nan(length(cdat.mean.f1), length(vowspprednames));
    allpred(:,1) = ones(length(cdat.mean.f1), 1);
    allpred(:,2) = zscore(cdat.mean.f1);
    allpred(:,3) = zscore(cdat.mean.f2);
    allpred(:,4) = allpred(:,2).*allpred(:,3);
    allpred(:,5) = cdat.mean.inVowSpace;
    allpred(:,6) = allpred(:,5).*allpred(:,2);
    allpred(:,7) = allpred(:,5).*allpred(:,3);
    allpred(:,8) = allpred(:,5).*allpred(:,4);
    
    %% run model
    for cm = 1:length(linmDefine2)
        cpred = allpred(:,linmDefine2{cm,2});
        curmfit.bnames = vowspprednames;
        curmfit.b = nan(size(y,1), size(cpred,2),size(y,2));
        curmfit.rsq = nan(cdat.elcount,1);
        curmfit.RSS = nan(cdat.elcount,1);
        curmfit.y_pred = nan(size(cy));
        curmfit.modelName = linmDefine2{cm,1};
        % cross-validated parameters
        [curmfit.y_pred(incltr,:),curmfit.b(:,:,incltr), curmfit.rsq, curmfit.RSS,curmfit.rCorr,curmfit.Ftest, curmfit.redM] =...
            crossvalFit_v3_vect(cpred(incltr,:),cy(incltr,:), 'linear');
        % write it out
        linmfits(csubj,cmnum+cm) = curmfit;
        clear curmfit;
    end
end

%% add info to models
% add subject and electrode to models
for csubj = 1:size(linmfits,1)
    for cm = 1:size(linmfits,2)
        linmfits(csubj,cm).subj = csubj*ones(size(linmfits(csubj,cm).rsq));
        linmfits(csubj,cm).el= (1:length(linmfits(csubj,cm).rsq));        
    end
end
% put F values and reduced model data in flat hierarchy in linmfits
for cm = 1:numel(linmfits)
    fns = fieldnames(linmfits(cm).Ftest);
    for fi= 1:length(fns)
    linmfits(cm).(fns{fi}) = linmfits(cm).Ftest.(fns{fi});
    end
    
    fns = fieldnames(linmfits(cm).redM);
    for fi= 1:length(fns)
    linmfits(cm).(fns{fi}) = linmfits(cm).redM.(fns{fi});
    end
end
%% concatenate info from models
linmAll.modnames = {linmfits(1,:).modelName}; 
linmAll.subj = [linmfits(:,1).subj];
linmAll.el = [linmfits(:,1).el];
linmAll.rsq=[];for i = 1:size(linmfits,2), linmAll.rsq(:,i) = [linmfits(:,i).rsq]; end
[linmAll.maxrsq,linmAll.maxmod] = max(linmAll.rsq,[],2);
% concatenate betas and unique variance for each model
for cm = 1:size(linmfits,2)
    linmAll.betas{cm}=[];
    linmAll.betaSD{cm}=[];
    for i = 1:size(linmfits,1)
        linmAll.betas{cm} = cat(1,linmAll.betas{cm},nanmean(linmfits(i,cm).b,3));         
        linmAll.betaSD{cm} = cat(1,linmAll.betaSD{cm},nanstd(linmfits(i,cm).b,0,3));         
    end    
    linmAll.betaName{cm} = linmfits(1,cm).bnames;
    linmAll.redMrsq{cm} = [linmfits(:,cm).uvar];
end
% unique variance
linmAll.uvnames = {'mainEff', 'F1', 'F2', 'IA'};%, 'F1/F2', 'F2/F1','IA & F1/F2 & F2/F1'};
linmAll.uv(:,1) = linmAll.rsq(:,1);
linmAll.uv(:,2) = linmAll.rsq(:,1) - linmAll.rsq(:,strcmpi(linmAll.modnames, 'F2'));
linmAll.uv(:,3) = linmAll.rsq(:,1) - linmAll.rsq(:,strcmpi(linmAll.modnames, 'F1'));
linmAll.uv(:,4) = linmAll.rsq(:,2) - linmAll.rsq(:,1); %linmAll.redMrsq{15}(4,:)';%

% unique variance for all models
for cm = 1:length(linmAll.modnames)
    linmAll.reduVar{cm} = [linmfits(:,cm).uvar]';
end
% F
for cm = 1:size(linmfits,2)
    linmAll.allF{cm} = [linmfits(:,cm).F];
end



%% plot predictor correlations
figure, imagesc(corr(allpred)), caxis([-1 1]); colormap(rbcm); colorbar;
xticks(1:7); xticklabels(allprednames);
yticks(1:7); yticklabels(allprednames);
num2str(corr(allpred),2)




%% --------- refit linear models on subset of stimuli within standard vowel space
for csubj = 1:length(alldat)
    cdat = alldat(csubj);
    %% hga - y
    y = squeeze(nanmean(cdat.mean.resp(:,cdat.befaft(1)+(10:40),:),2));
    cy = y';
%     cy(isfinite(cy(:,1)),:) = zscore(cy(isfinite(cy(:,1)),:),0,1);
    incltr = isfinite(cy(:,1)) & cdat.mean.inVowSpace;
    cy(incltr,:) = zscore(cy(incltr,:),0,1);
    %% all predictors
    allpred = nan(sum(incltr), length(allprednames));
    allpred(:,1) = ones(length(cdat.mean.f1(incltr)), 1);
    allpred(:,2) = zscore(cdat.mean.f1(incltr));
    allpred(:,3) = zscore(cdat.mean.f2(incltr));
    allpred(:,4) = allpred(:,2).*allpred(:,3);
    allpred(:,5) = zscore(cdat.mean.f1(incltr)./cdat.mean.f2(incltr));
    allpred(:,6) = zscore(cdat.mean.f2(incltr)./cdat.mean.f1(incltr));
    allpred(:,7) = zscore(cdat.mean.f2(incltr)-cdat.mean.f1(incltr));
    
    %% all models, predictor selection and model fit
    for cm  = 1:size(linmDefine,1)
        cpred = allpred(:,linmDefine{cm,2});
        curmfit.bnames = allprednames(linmDefine{cm,2});
        curmfit.b = nan(size(y,1), size(cpred,2),size(y,2));
        curmfit.rsq = nan(cdat.elcount,1);
        curmfit.RSS = nan(cdat.elcount,1);
        curmfit.y_pred = nan(size(cy));
        curmfit.modelName = linmDefine{cm,1};
        % cross-validated parameters
        [curmfit.y_pred(incltr,:),curmfit.b(:,:,incltr), curmfit.rsq, curmfit.RSS,curmfit.rCorr,curmfit.Ftest, curmfit.redM] =...
            crossvalFit_v3_vect(cpred,cy(incltr,:), 'linear');
        % write it out
        linmfitSm(csubj,cm) = curmfit;
        clear curmfit;
    end
end
% add subject and electrode to models
for csubj = 1:size(linmfitSm,1)
    for cm = 1:size(linmfitSm,2)
        linmfitSm(csubj,cm).subj = csubj*ones(size(linmfitSm(csubj,cm).rsq));
        linmfitSm(csubj,cm).el= (1:length(linmfitSm(csubj,cm).rsq));        
    end
end
% put F values and reduced model data in flat hierarchy in linmfitSm
for cm = 1:numel(linmfitSm)
    fns = fieldnames(linmfitSm(cm).Ftest);
    for fi= 1:length(fns)
    linmfitSm(cm).(fns{fi}) = linmfitSm(cm).Ftest.(fns{fi});
    end
    
    fns = fieldnames(linmfitSm(cm).redM);
    for fi= 1:length(fns)
    linmfitSm(cm).(fns{fi}) = linmfitSm(cm).redM.(fns{fi});
    end
end
% concatenate info from models
linmAllSm.modnames = linmDefine(:,1);
linmAllSm.subj = [linmfitSm(:,1).subj];
linmAllSm.el = [linmfitSm(:,1).el];
linmAllSm.rsq=[];for i = 1:size(linmfitSm,2), linmAllSm.rsq(:,i) = [linmfitSm(:,i).rsq]; end
[linmAllSm.maxrsq,linmAllSm.maxmod] = max(linmAllSm.rsq,[],2);

% concatenate betas and unique variance for each model
for cm = 1:size(linmfitSm,2)
    linmAllSm.betas{cm}=[];
    for i = 1:size(linmfitSm,1)
        linmAllSm.betas{cm} = cat(1,linmAllSm.betas{cm},nanmean(linmfitSm(i,cm).b,3));         
    end    
    linmAllSm.betaName{cm} = linmfitSm(1,cm).bnames;
    linmAllSm.redMrsq{cm} = [linmfitSm(:,cm).uvar];
end
%% unique variance

linmAllSm.uvnames = {'mainEff', 'F1', 'F2', 'IA'};%, 'F1/F2', 'F2/F1','IA & F1/F2 & F2/F1'};
linmAllSm.uv(:,1) = linmAllSm.rsq(:,1);
linmAllSm.uv(:,2) = linmAllSm.rsq(:,1) - linmAllSm.rsq(:,strcmpi(linmAllSm.modnames, 'F2'));
linmAllSm.uv(:,3) = linmAllSm.rsq(:,1) - linmAllSm.rsq(:,strcmpi(linmAllSm.modnames, 'F1'));
linmAllSm.uv(:,4) = linmAllSm.rsq(:,2) - linmAllSm.rsq(:,1); 


% F
for cm = 1:size(linmfitSm,2)
    linmAllSm.allF{cm} = [linmfitSm(:,cm).F];
end

