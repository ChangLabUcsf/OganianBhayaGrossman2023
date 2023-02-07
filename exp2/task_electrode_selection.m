% all electrodes in analysis
AnaElTask = max(linmAll.rsq(:,1:2),[],2) >.1;
sum(AnaElTask)

%% excl EC219, patient with epilepsy focus in speech cortex
% AnaElTask = AnaElTask & linmAll.subj'~=2;

%% for neuron revision 
% repeating all analyses in here with electrodes from spanish monolinguals
% only EC214, EC225
% spanishSid = [1 5];
% AnaElTask = AnaElTask & ismember(linmAll.subj', spanishSid); 
% sum(AnaElTask)
%% electrode color by F1/F2 beta
elcols=[];
colCoord = min(ceil((ecog_norm(linmAll.betas{2}(AnaElTask, [2 3]), 1)+1.0001)*100),200);
for i = 1:sum(AnaElTask), elcols(i,:) = cm2d(colCoord(i,1),colCoord(i,2),:);end
