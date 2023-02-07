
% "STARTUP.M" NEEDS TO BE RUN BEFORE THE FOLLOWING CODE

%% ------------ Vowel Formant Encoding in Natural Speech ------------------

[inflections]=loadInfl(betaInfo, SIDs, 1, Dvow, bef, aft);

pitchSel = ''; % '', 'gt170', 'lt170'

decode_model = ['decode_vowel_Spanish_20-50_allVowel_vowelftestEl-single_r0.15_' pitchSel '.mat'];
load([datapath '/ecog_decode/' decode_model])

neur_model{1} = neur;
acs_model{1} = acs;
          
sids = neur.respidx(:, 1);
els = neur.respidx(:, 2);

ssids = cellfun(@(x) str2double(x(3:end)), SIDs);
idx = ismember(sids, ssids);

neur_model{1}.decode_single = table(sids(idx), els(idx), neur.confMat(idx)', ...
    'VariableNames',  {'SID', 'el', 'conf'});

% for pooled within subj: 'decode_vowel_Spanish_20-50_allVowel_vowelftestEl-subj_r0.15_'

% only 5 subjects used to maximize overlapping trials
% for pooled across subjs: 'decode_vowel_Spanish-mini_20-50_allVowel_vowelftestEl_r0.15_'
% with EC214: 'decode_vowel_Spanish-mini3_20-50_allVowel_vowelftestEl_r0.15_'
% without EC214: 'decode_vowel_Spanish-mini2_20-50_allVowel_vowelftestEl_r0.15_'
decode_model = ['decode_vowel_Spanish-mini2_20-50_allVowel_vowelftestEl_r0.15_' pitchSel '.mat'];
load([datapath '/ecog_decode/' decode_model])

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* inflections ...
    formant *model 

% SVM: Run decoding on sub-populations

% bootstrap reps
reps = 50;
neuresp = cell2mat(inflections.resp);
tps = bef:bef+30;

% cases where both formants are encoded significantly in rsq
joint_encoding = inflections.sr_f1>0 & inflections.sr_f2>0;

% decoding single

% for each electrode, confusion matrix is 5x5xtps 
decode_single = neur_model{1}.decode_single;
confSingleNaN = arrayfun(@(x) max(decode_single.conf{x}, [], 3, 'omitnan'), ...
    1:height(decode_single), 'UniformOutput',false);
confSingleNaN = cat(3, confSingleNaN{:});

% find best single electrode per vowel category
selel = nan(5, 1);
for v = 1:5
    % all pairwise accuracies including this pair
    tmp = [squeeze(confSingleNaN(v, :, :)); squeeze(confSingleNaN(:, v, :))];
    [~, els] = maxk(mean(tmp, 'omitnan'), 15);
    els = squeeze(els);

    ctr = 1;
    while isnan(selel(v))
        % find an electrode in inflections and decode single
        SID = decode_single.SID(els(ctr));    
        el = decode_single.el(els(ctr));   

        % find overlap in inflections and decode_single
        idx = ismember(inflections.SID, SID) & ismember(inflections.el, el);
        if sum(idx) && ~ismember(find(idx), selel)
            selel(v) = find(idx);
        end
        ctr = ctr+1;
    end
    clear ctr
end

confMat = cell(reps, 5);
pconfMat = cell(reps, 5);
tpsfrmOnset = nan(5, 5);
   
incl_vow = {'a', 'e', 'i', 'o', 'u'};
label = {'F1-/F2+', 'F1+/F2-', 'Joint'};
betas = {[], [], [], [], []};

filename_single = ...
    'ecog_decode/decode_vowel_Spanish_maxFstattp_allVowel_vowelftestEl-single_.mat';
filename_pop = ...
    'ecog_decode/decode_vowel_Spanish_maxFstattp-pca_allVowel_vowelftestEl-both.mat';
% filename_pop = ... % all
%     'ecog_decode/decode_vowel_Spanish_maxFstattp-pca_allVowel_vowelftestEl_.mat';

% if the decoding has not yet been run
if ~exist([datapath filename_single], 'file') 

    figure; 
    tic
    for vowel = 1:5   % only run electrodes with best mean performance for this vowel 
        
        X = inflections.resp{selel(vowel)}(1, tps, :);
        [confMat(:, vowel), pconfMat(:, vowel), beta, wind, nanel, ~, classes] = ...
            decode_permute(X, Dvow.vowelType, Dvow.vowel, [], reps, 0); 
        betas(vowel) = {beta};
        tpsfrmOnset(vowel, :) = wind;
    
        % visualization
        subplot(1, 5, vowel);
        for i = 1:10 
            r = rand(reps, 1);
            a = cell2mat(confMat(:, vowel)); 
            hold on; 
            scatter(r/4+i, a(:, i), 55, 'filled'); 
            b = cell2mat(pconfMat(:, vowel)); 
            scatter(r/4+i, b(:, i), 55, [0.6 0.6 0.6], 'filled', ...
                'MarkerFaceAlpha', 0.3);           
        end
        xticks(1:10);
        yline(0.5);
        xticklabels(arrayfun(@(x) strjoin(classes{x}, '-'), 1:10, ...
            'UniformOutput', false));
        title(incl_vow{vowel});
    end
    toc
    save([datapath filename_single], ...
        'selel', 'confMat', 'pconfMat', 'betas',  'tpsfrmOnset');
    clear selec confMat pconfMat betas tpsfrmOnset
end
    
if ~exist([datapath filename_pop], 'file')
    % decoding subpopulations
    confMat = cell(reps, 3);
    pconfMat = cell(reps, 3);
    
    % first in  f1-, f1+,, all joint
    f1_minus = find(joint_encoding & inflections.sl_f1<1 & inflections.sl_f2>0);
    f1_plus = find(joint_encoding & inflections.sl_f1>0 & inflections.sl_f2<1);
    elecs = {f1_minus, f1_plus, union(f1_minus,  f1_plus)};
    
    betas = {[], [], []};
    selec = {[], [], []};
    tpsfrmOnset = nan(3, 5);
    for cond = 1:3 % conditions: f1+, f1-, all joint
        tic
        X = neuresp(elecs{cond}, tps, :);           
        disp(['Number of electrodes: ' num2str(length(elecs{cond}))])
    
        [confMat(:, cond), pconfMat(:, cond), beta, wind, nanel, ...
            ~, classes] = decode_permute(X, Dvow.vowelType, Dvow.vowel, ...
            [], reps, 1);% 0.25
        betas(cond) = {beta};
        % all time points used per population classifier
        tpsfrmOnset(cond, :) = wind;
        selec{cond} = elecs{cond}(~nanel);
    
        % visualization
        figure; 
        for i = 1:10 
            r = rand(reps, 1);
            a = cell2mat(confMat(:, cond)); 
            hold on; 
            scatter(r/4+i, a(:, i), 55, 'filled'); 
            b = cell2mat(pconfMat(:, cond)); 
            scatter(r/4+i, b(:, i), 55, [0.6 0.6 0.6], 'filled', ...
                'MarkerFaceAlpha', 0.3);           
        end
        xticks(1:10);
        yline(0.5);
        xticklabels(arrayfun(@(x) strjoin(classes{x}, '-'), 1:10, ...
            'UniformOutput', false));
        title(label{cond})
        toc
    end
    
    % save out results
    save([datapath filename_pop], ...
        'selec', 'confMat', 'pconfMat', 'betas', 'tpsfrmOnset');
    clear selec confMat pconfMat betas tpsfrmOnset
end

all_decode = load([datapath filename_pop]);
single_decode = load([datapath filename_single]);

% aggregate all useful information for plotting
% reps by electrode condition
acc = nan(10, reps, 4);
perm_acc = nan(10, reps, 4);
load('svmClassOrder.mat');

% pair by pair by electrode condition (only single electrode and all)
confAll = nan(5, 5, 2, reps);
ctr = 1;
for v1 = 0:4
    for v2 = v1+1:4
        
        % find index of pair
        pair = find(all([classOrd.classesNum] == [v1; v2]));
        
        % scatter/boxplot single electrode best accuracy
        tmp_acc = cellfun(@(x) x(pair), single_decode.confMat);        
        [~, elec] = max(mean(tmp_acc, 1));
        acc(ctr, :, 1) = tmp_acc(:, elec)';
        confAll(v1+1, v2+1, 1, :) = tmp_acc(:, elec)';  
        confAll(v2+1, v1+1, 1, :) = tmp_acc(:, elec)';  
        
        tmp_pacc = cellfun(@(x) x(pair), single_decode.pconfMat);        
        perm_acc(ctr, :, 1) = tmp_pacc(:, elec)';
        clear tmp*
        
        % scatter/boxplot subsets/all with accuracy
        acc(ctr, :, 2:end) = cellfun(@(x) x(pair), all_decode.confMat);                
        perm_acc(ctr, :, 2:end) = cellfun(@(x) x(pair), all_decode.pconfMat); 

        % use all electrodes (~56) versus the miniset
        confAll(v1+1, v2+1, 2, :) = acc(ctr, :, 4);
        confAll(v2+1, v1+1, 2, :) = acc(ctr, :, 4);
        ctr = ctr+1;
    end   
end

decode_all.confAll = confAll;
decode_all.acc = acc;
decode_all.perm_acc = perm_acc;
decode_all.class_ord = classOrd;
decode_all.condition = label;
decode_all.reps = reps;
decode_all.single = single_decode;

% copy over
decode_all.confMat = all_decode.confMat;
decode_all.selec = all_decode.selec;
decode_all.betas = all_decode.betas;
decode_all.tpsfrmOnset = all_decode.tpsfrmOnset;
clear all_decode

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections* formant *model

%% Figure 4 (a 1) - averaged pairwise comparisons

confAll = decode_all.confAll;
colors = getColors(1);
reps = decode_all.reps;

figure;
pairavg = nan(6, 50);
pairavg(6, :) = mean(cell2mat(decode_all.confMat(:, 3)), 2);
for v = 1:5
    % average across all pairs for every rep
    pairavg(v, :) = mean(cell2mat(decode_all.single.confMat(:, v)), 2);
    boxchart(ones(50, 1)*v, squeeze(pairavg(v, :))', 'MarkerColor', 'none', 'LineWidth', 0.5, ...
        'BoxFaceColor',colors(v, :)); hold on;

    [~, p] = ttest2(pairavg(6, :), pairavg(v, :));
    med_diff = (median(pairavg(6, :))-median(pairavg(v, :)))*100;
    disp(['Median difference: ' num2str(med_diff) ', pval: ' num2str(p)]);
end

boxchart(ones(50, 1)*7, pairavg(6, :)', 'MarkerColor', 'none', ...
    'LineWidth', 0.5, 'BoxFaceColor',[0.1 0.1 0.1]); hold on;

% permutation baseline
pacc = [cellfun(@(x) mean(x), decode_all.single.pconfMat),...
    mean(decode_all.perm_acc(:, :, 4), 1)'];
pacc_tmp = [pacc(:, 1), pacc,  pacc(:, 6)];
shadedErrorBar(0.5:7.5, mean(pacc_tmp), std(pacc_tmp), ...
    'patchSaturation', 0.02);


set(gca, 'FontSize', 15);
xlim([0 8]);
xticks([1:5 7])
xticklabels({'E1', 'E2', 'E3', 'E4', 'E5', 'all'});
yticks([0.5, 0.7]);
yticklabels({'50', '70'});
ylabel('% correct across pairs');

% stats on single versus population
confAllNan = confAll;
for p = 1:reps
    for r = 1:5
        for c = r+1:5
            confAllNan(c, r, :, p) = NaN;
        end
    end
end
pairs = reshape(confAllNan, [5*5, 2, reps]);
pairs(isnan(pairs(:, 1, 1)), :, :) = [];

p = nan(10, 1);
for i = 1:10
    [~, p(i)] = ttest2(squeeze(pairs(i, 1, :))', squeeze(pairs(i, 2, :))');
end
disp(['Pairs at p<0.0001 = ' num2str(sum(p<0.0001)) '/10']);
[p, ~, stats] = ranksum(reshape((mean(pairs(:, 2, :), 3)), [], 1), ...
    reshape((mean(pairs(:, 1, :), 3)), [], 1));
disp(['All pair p = ' num2str(p) ' ']);

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections* formant *model


%% Figure 4 (a 2) - ERPs for selected single electrode decoding 
fig = figure('Position',[0 100 200 800]);
incl_vow = {'a', 'e', 'i', 'o', 'u'};
els = decode_all.single.selel; % 
for ctr = 1:5  

    SID = ['EC' num2str(inflections.SID(els(ctr)))];
    elec = inflections.el(els(ctr));

    Dvow = addtoDD(Dvow, 'dimex', bef, aft, {SID});
   
    subplot(5, 1, ctr)

    plotVowelErp(Dvow, SID, elec, fig, getColors(1), bef./100);
    ylabel('HGA (norm)');
    xlabel('Time (s)');
    set(gca, 'FontSize', 15);
    yticks(-0.5:0.5:1); 
    title(incl_vow(ctr));
end

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections* formant *model

%% Figure 4 (c) - pairwise decoding comparisons

load('svmClassOrder.mat');
ord = [classOrd.classesNum];

% summary with accuracy
acc = decode_all.acc;
pacc = decode_all.perm_acc;

figure;
% pair by pair by electrode condition (only single electrode and all)
ctr = 1;
for v1 = 0:4
    for v2 = v1+1:4
        subplot(4, 4, gridCount(v1+1, v2, 4, 4));
        
        % find index of pair
        pair = find(all([decode_all.class_ord.classesNum] == [v1; v2]));
        
        boxchart(squeeze(acc(pair, :, :)), 'MarkerColor', 'none', 'LineWidth', 0.5); 
        hold on;
        
        pacc_tmp = [squeeze(pacc(pair, :, 1))' squeeze(pacc(pair, :, 1:4)) ...
            squeeze(pacc(pair, :, 4))'];
        shadedErrorBar(0:5, mean(pacc_tmp), std(pacc_tmp), ...
            'patchSaturation', 0.04);
        
        plot(mean(squeeze(acc(ctr, :, :)), 'omitnan'), '-k', 'LineWidth', 1.25);
        xticklabels({'s', '-/+', '+/-', 'b'});
        
        % check for significant difference between two subpopulations
        % F1-/F2+ and F1+/F2-
        if size(acc, 3)>3
            [~, p, ~] = ttest2(acc(pair, :, 2), acc(pair, :, 3));
        else
            [~, p, ~] = ttest2(acc(pair, :, 1), acc(pair, :, 2));
        end
        if p<0.0001 
            plot([2, 3], [0.75 0.75], '-k');
            text(2.3, 0.76, 'p<.0001', 'FontSize', 12);
        end  
        
        ylim([0.43 0.81]);
        xlim(categorical([1 size(acc, 3)]));
        if ismember(gridCount(v1+1, v2, 4, 4), [1, 6, 11, 16])
            yticks([0.5, 0.80]);
            yticklabels({'50', '80'});
        else
            yticks([]);
        end

        if gridCount(v1+1, v2, 4, 4)~=16
            xticks([]);
        end
        title(strjoin(decode_all.class_ord(pair).classes)) 
        ctr = ctr+1;
    end
end

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections* formant *model

%% Figure 4 (d) - plot RFs of two electrode subsets

figure; 
ax(1) = subplot(2, 2, 1);
elidx = inflections(decode_all.selec{1}, :);
elidx.SID = arrayfun(@(x) ['EC' num2str(x)], elidx.SID, ...
    'UniformOutput', false);

vals = plotAvgImagescResponse(Dvow, elidx, getColors(1), [], ...
    rescale(max([elidx.sr_f1, elidx.sr_f2], [], 2)));
title('F1-/F2+');
set(gca, 'FontSize', 15);
caxis([-0.05 0.65])
% cmap = brewermap(256, 'PuOr'); 
% colormap([1 1 1; flipud(cmap(40:240, :))]);
box off;
yticks([]);
xticks([]);
xlim([800, 3000]);
ylim([200 800]);
brighten(0.2);

subplot(2, 2, 2);
scatter(elidx.sr_f1, elidx.sr_f2, 55, 'k', 'filled');
h=refline(1, 0);
yticks([0 1]);
xticks([0 1]);
xlabel('F1 R^2');
ylabel('F2 R^2');
set(gca, 'FontSize', 15);
h.Color = 'k';
h.LineStyle = '--';

ax(2) = subplot(2, 2, 3);
elidx = inflections(decode_all.selec{2}, :);
elidx.SID = arrayfun(@(x) ['EC' num2str(x)], elidx.SID, ...
    'UniformOutput', false);

vals2 = plotAvgImagescResponse(Dvow, elidx, getColors(1), [], ...
    rescale(max([elidx.sr_f1, elidx.sr_f2], [], 2)));
title('F1+/F2-');
set(gca, 'FontSize', 15);
caxis([-0.05 0.65]);
box off;
yticks([]);
xticks([]);
xlim([800, 3000]);
ylim([200 800]);
brighten(0.2);

subplot(2, 2, 4);
scatter(elidx.sr_f1, elidx.sr_f2, 55, 'k', 'filled');
h=refline(1, 0);
yticks([0 1]);
xticks([0 1]);
xlabel('F1 R^2');
ylabel('F2 R^2');
set(gca, 'FontSize', 15);
h.Color = 'k';
h.LineStyle = '--';

% marking the difference lines on the graphs
medoids = plotVariance(Dvow.formantVals, Dvow.vowel, [], 0, 0);
diffacc = nan(5, 5);
for v1 = 1:5
    for v2 = v1+1:5
        tmp = [abs(triu(vals-vals')) abs(triu(vals2-vals2'))];
        tmp = tmp(tmp>0);

        subplot(2, 2, 1)
        width = discretize(abs(diff(vals([v1, v2]))), ...
            linspace(min(tmp), max(tmp), 8))*0.6;
        plot([medoids(v1, 2) medoids(v2, 2)], ...
            [medoids(v1, 1) medoids(v2, 1)], 'Color', 'k', ...
            'LineWidth', width);

        subplot(2, 2, 3)
        width = discretize(abs(diff(vals2([v1, v2]))), ...
            linspace(min(tmp), max(tmp), 8))*0.6;
        plot([medoids(v1, 2) medoids(v2, 2)], ...
            [medoids(v1, 1) medoids(v2, 1)], 'Color', 'k', ...
            'LineWidth', width);
        diffacc(v1, v2) = abs(diff(vals([v1, v2]))) - abs(diff(vals2([v1, v2])));
    end
end

% figure;
% tmp = diffacc;
% tmp(isnan(tmp)) = 0;
% imagesc(tmp);

% formatting
% xticks(1:5); yticks(1:5);
% xticklabels({'a', 'e', 'i', 'o', 'u'});
% yticklabels({'a', 'e', 'i', 'o', 'u'});
% colormap(flipud(brewermap(15, 'RdBu')));
% caxis([-0.1 0.1]);
% ylim([0.5 4.5]);
% xlim([1.5 5.5]);

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections* formant *model
%% Figure 4 (e) - summary of subpopulation decoding

load('svmClassOrder.mat');
ord = [classOrd.classesNum];

% find index of pair
% summary with accuracy
acc = arrayfun(@(i) cellfun(@(x) x(all([classOrd.classesNum] == ord(:, i))), ...
    decode_all.confMat), 1:size(ord, 2), 'UniformOutput', false);  
acc = cat(3, acc{:});

summaryDecode(Dvow, acc, 0.0001)

confacc = nan(3, 5, 5);
ctr = 1;
for pair = ord
    confacc(:, pair(1)+1, pair(2)+1) = mean(acc(:, :, ctr));
    ctr = ctr+1;
end

figure;
tmp = squeeze(confacc(1, :, :)-confacc(2, :, :));
tmp(isnan(tmp)) = 0;
imagesc(tmp);

% formatting
xticks(1:5); yticks(1:5);
xticklabels({'a', 'e', 'i', 'o', 'u'});
yticklabels({'a', 'e', 'i', 'o', 'u'});
colormap(flipud(brewermap(15, 'RdBu')));
caxis([-0.1 0.1]);
ylim([0.5 4.5]);
xlim([1.5 5.5]);


clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections* formant *model

%% Figure 4 - SUPP: STATS on subpopulation differences

% each row is a pair, first column is median accuracy, second column is t-stat, 
% third is p val
supp_stats =  array2table(zeros(0,4), 'VariableNames', ...
    {'pair', 'acc', 'z-val', 'p-val'});
for i = 1:size(decode_all.acc)
    Y = [squeeze(decode_all.acc(i, :, 2)); squeeze(decode_all.acc(i, :, 3))];
    %[p, h, stats] = ranksum(Y(1, :), Y(2, :));
    [h,p, ~, stats] = ttest2(Y(1, :), Y(2, :));
    pair_label = strjoin(decode_all.class_ord(i).classes);
    supp_stats = [supp_stats; {pair_label, mean(Y, 2), stats.tstat, p}];
end
disp(supp_stats)

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections* formant *model supp_stats
%% Figure 4 - OPT: logistic regression code

resp = cell2mat([inflections.resp]);
AUC = nan(3, 5, 5, 50);

% cases where both formants are encoded significantly in rsq
% joint_encoding = inflections.sr_f1>0 & inflections.sr_f2>0;
% f1_minus = find(joint_encoding & inflections.sl_f1<1 & inflections.sl_f2>0);
% f1_plus = find(joint_encoding & inflections.sl_f1>0 & inflections.sl_f2<1);
% use same electrodes as in the svm case
elecs = decode_all.selec;

for v = 1:4
    for v2 = v+1:5
        for type = 1:3 % F1-, F1+, both
            % all vowels in these categories
            idx = Dvow.vowelType == v | Dvow.vowelType == v2;
            tps = 35:40; % close to time points used for svm classifiers
            X = reshape(resp(elecs{type}, tps, idx), ...
                length(elecs{type})*length(tps), []);
        
            % remove columns and rows where data is missing
            % trials for which over electrodes do not have data
            % remanining electrodes with missing trials
            % A is electrode/time x trials
            nancol = sum(isnan(X))>0.1*size(X, 1);
            X(:, nancol) = [];
            nanrow = sum(isnan(X), 2)>0;
            X(nanrow, :) = [];
            assert(sum(nanrow)==0)
        
            y = double(Dvow.vowelType(idx)==v);
            y(nancol) = [];
    
            [~, score, ~, ~, exp] = pca(X'); 
            n_comp = find(diff(cumsum(exp)>90));           

            for n = 1:50
                % sample 200 of each class
                c1 = find(y==1); c2 = find(y==0);
                idx = [c1(randperm(length(c1), 110)), c2(randperm(length(c2), 110))];

                % perform logistic regression on each subset and all electrodes
                CVMdl = fitclinear(score(idx, 1:n_comp), y(idx),'KFold',5,...
                    'Learner','logistic','Regularization','lasso');
                labels = kfoldPredict(CVMdl);
                AUC(type, v, v2, n) = sum(labels==y(idx)')/length(idx);          
            end
        end
    end
end

AUC_all = AUC;

% logistic regression visualization
load('svmClassOrder.mat');
ord = [classOrd.classesNum];

% find index of pair
figure;
a = squeeze(mean(AUC_all, 4, 'omitnan'));
tmp = squeeze(a(1, :, :)-a(2, :, :));
tmp(isnan(tmp))=0;
imagesc(tmp);

% formatting
xticks(1:5); yticks(1:5);
xticklabels({'a', 'e', 'i', 'o', 'u'});
yticklabels({'a', 'e', 'i', 'o', 'u'});
colormap(flipud(brewermap(15, 'RdBu')));
caxis([-0.1 0.1]);
ylim([0.5 4.5]);
xlim([1.5 5.5]);

% plot summary
% summary with accuracy
acc_auc = arrayfun(@(i) squeeze(AUC_all(:, ord(1, i)+1, ord(2, i)+1, :)), ...
    1:size(ord, 2),'UniformOutput', false); 
acc_auc = permute(cat(3, acc_auc{:}), [2, 1, 3]);
summaryDecode(Dvow, acc_auc, 0.0001);


%% Figure 4 (a, b) - Decoding from single and across electrodes

% for all single electrodes
decode_single = neur_model{1}.decode_single;
confSingle = nan(5, 5, height(decode_single));
weights = nan(height(decode_single), 5);
    
% for each electrode
for i = 1:height(decode_single)
    conf = decode_single.conf{i};
    %[max_acc(i), tps(i)] = max(squeeze(nanmean(max(conf), 2)));
    c = max(conf, [], 3, 'omitnan');   
    c(isnan(c)) = 0;
    confSingle(:, :, i) = c+c';
    
    % find vowel selectivity index - average pairwise acc for single vowel
    c = confSingle(:, :, i);
    c(c==0) = NaN; 
    weights(i, :) = mean(c, 'omitnan');
end

% add weights to table
decode_single.weights = weights;
neur{1}.decode_single = decode_single;

% for all electrodes together
confAll = decode_all.confMat{1};
%[max_acc(i+1), tps(i+1)] = max(squeeze(nanmean(max(confAll), 2)));
%range = max(1, tps(i+1)-wind):min(size(conf, 3), tps(i+1)+wind);

% average single electrode confusion matrices vs. all electrode
% previous case, where look at best single electrode accuracy
% confs = {max(confSingle, [], 3, 'omitnan'), max(confAll, [], 3, 'omitnan')};

% new case, where pick best electrode per row
confMax = nan(5, 5);
el = nan(5, 1);
for v = 1:5
    % maximum electrode per row
    confSingle(confSingle==0)= NaN;
    [~, el(v)] = max(mean(confSingle(v, :, :), 'omitnan'), [], 3, 'omitnan'); 
    confMax(v, :) = confSingle(v, :, el(v));
end
% selected electrodes per category
selel = decode_single(el, :);
confs = {confMax, max(confAll, [], 3, 'omitnan'), ...
    max(acs_model{1}.confMat{1}, [], 3, 'omitnan')};
clear range conf

% confusion matrix
figure;
for t = 1:2
    
    incl_vows = neur_model{1}.incl_vow;
    confMat = confs{t};
    confMat(isnan(confMat))=0;
    if t == 2
        subplot(5, 2, [2 4 6 8 10]);
        imagesc(confMat+confMat');
        caxis([0.5 1]);
                
        colormap(gca, brewermap(30, 'YlGn'));
        hcb=colorbar;
        hcb.Ticks=0.5:0.1:1;
        ylabel(hcb, 'pairwise classification accuracy');
        
        xticks(1:length(incl_vows)); yticks(1:length(incl_vows));
        xticklabels(incl_vows);
        yticklabels(incl_vows);
        
        for x = 1:5
            for y = x+1:5
                text(y-0.3, x, num2str(confMat(x, y), '%.2f'), ...
                    'FontSize', 15);
            end
        end
    else
        for v = 1:5
            subplot(5, 2, 1+(v-1)*2);       
            imagesc(confMat(v, :));
            caxis([0.5 1]);
            yticks(1);
            yticklabels(incl_vows(v))
            
            colormap(gca, brewermap(30, 'YlGn'));
            for v2 = find(1:5 ~= v)
                text(v2-0.3, 1, num2str(confMat(v, v2), '%.2f'), ...
                    'FontSize', 15);
            end
            
            if v == 5
                xticklabels(incl_vows);
            else
                xticks([]);
            end
            set(gca,'fontsize',15);
        end
    end       
    set(gca,'fontsize',15);       
end

% revisualize as two clouds
figure;
subplot(1, 3, [1 2]);
cols = getColors(1);
for v = 1:5    
    for i = 1:2
        c = confs{i};
        if i==2
            c(isnan(c))=0;
            c = c+c';            
        end
        c(c==0) = NaN;
        y(i) = mean(c(v, :), 'omitnan');
    end
    plot(y, 'Color', cols(v, :)); hold on;
    scatter([1, 2], y, 85, cols(v, :), 'filled', 'MarkerFaceAlpha', 0.5);  
end

xlim([0.5 2.5]);
ylim([0.55 0.85]);
xticks([1 2]);
yticks(0.5:0.1:1);
xticklabels({'single', 'pooled'});
xlabel('Decoding Type');
ylabel('Accuracy')
set(gca, 'FontSize', 15);

subplot(1, 3, 3);
boxplot(confs{2}(:) - confs{1}(:));
xticks([]);
title('Pooled - Single');
box off;

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections* formant *model el
%% Figure 4 (b) - ERPs for population level decoding (run cell before)

decode_single = neur_model{1}.decode_single;

numvows = 5;
maxtp = nan(5, 1);
% find maximum time point per svm
svmOut = decode_all.svmOut;
classes = [svmOut.classesNum];
classes = mat2cell(classes, 2, ones(10, 1));
% comps = find(sum(classes==v-1)>0);

% length(tps)
allresp = nan(size(decode_all.respidx, 1), bef+aft+1, numvows, numvows);
for el = 1:size(decode_all.respidx, 1)
    confmat = nan(numvows, numvows, length(tps));
        
    % ADDED: multiply by high gamma amplitude
    sid = decode_all.respidx(el, 1);
    elnum = decode_all.respidx(el, 2);
    
    % translate electrode index from betas to decode single table
    sel =  find(decode_single.SID == sid & decode_single.el == elnum);
      
    % hg_weights = nan(numvows, length(tps));
    hg_weights = nan(numvows, bef+aft+1);
    % add resp to decode_single
    for v = 1:numvows
        SID = ['EC' num2str(sid)];
        Dvow = addtoDD(Dvow, 'dimex', bef, aft, {SID});
        hg_weights(v, :) = mean(Dvow.(SID).resp(elnum, ...
            :, Dvow.vowelType==v), 3, 'omitnan'); %tps     
        
    end
    
    for v = 1:numvows
        allresp(el, :, v, v) = hg_weights(v, :)'.* ...
            max(squeeze(decode_single.weights(sel, v, :))); % one time point
            %squeeze(decode_single.weights(sel, v, :));
        for v2=1:5
            if v ~= v2
                svmOut = decode_all.svmOut;
                
                % find the comparison with weights
                c = cellfun(@(x) all(x==[v-1; v2-1] | x==[v2-1; v-1]), classes);
                if all(classes{c} == [v-1; v2-1])
                    w = svmOut(c).beta(:, el); 
                else
                    w = -1*svmOut(c).beta(:, el); 
                end

                [~, idx] = max(abs(w));
                w = w(idx);
                allresp(el, :, v, v2) = hg_weights(v2, :)'.*w.*20;
                clear c w
            end
        end
    end
    
    decode_single.hg_resp(sel) = {hg_weights};        
end

% find maximum svm weight per electrode
Dvow.all.resp = nan(1, length(tps), numvows);

figure('Renderer', 'painters', 'Position', [200 800 200 350*4]);
% plot weighted average ERPs
cols = getColors(1);
for v = 1:numvows
    subplot(5, 1, v);
    
    for v2 = 1:numvows
        % x = (tps-bef)./100;
        x = (-bef:aft)./100;
        if v2==v
            y = squeeze(mean(allresp(:, :, v, v2), 'omitnan'));
            plot(x, y-y(1), 'LineWidth', 2.5, 'Color', cols(v2, :)); hold on;
        else
            plot(x, squeeze(mean(allresp(:, :, v, v2), 'omitnan')), 'LineWidth', 2, ...
                'Color', cols(v2, :), 'LineStyle', ':'); hold on;
        end
    end
    xticks([]);
    yticks([]);
    ylim([-0.1 0.45]);
    xlim([-0.2 0.5]);
    xline(0);
    yline(0);
    box off;
end
yticks(0:0.2:0.4);
xticks([-0.2 0 0.4]);
xlabel('Time (s)');
set(gca, 'FontSize', 15);
box off;

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections* formant *model

%% Figure 4 - OPT: Decoding from single direction of encoding

%find direction of encoding
pf1nf2 = inflections(inflections.sr_f1>0 & inflections.sr_f2>0 & ...
     inflections.sl_f1>0 & inflections.sl_f2<1, :);
nf1pf2 = inflections(inflections.sr_f1>0 & inflections.sr_f2>0 & ...
    inflections.sl_f1<1 & inflections.sl_f2>0, :);
allsl = inflections(inflections.sr_f1>0 & inflections.sr_f2>0, :);

types = {pf1nf2, nf1pf2, allsl};
labels = {'f1+/f2-', 'f1-/f2+', 'all', 'single'};
data = Dvow;
ctr = 1;

neur = struct();
neur.tps = 20:50;
neur.frametype = 'all';
filename = [datapath ...
    '/ecog_decode/decode_vowel_Spanish_' num2str(neur.tps(1)) '-' num2str(neur.tps(end))...
    neur.frametype '_allVowel_slopeDirEl_r0.15_.mat'];

if ~isfile(filename)
    figure;
    confMat = cell(1, 3);
    svmOut = cell(1, 3);   
    
    neur.electype = '';
    neur.eleclabel = 'slopeDir';
    neur.incl_vow = {'a', 'e', 'i', 'o', 'u'};
    neur.subj = SIDs;
       
    neur.pca = 1;
    neur.respidx = {[], [], []};
    for t = types
        tbl = t{1};
        data.resp = [];
        respidx = [];
        for s = unique(tbl.SID)'
            % aggregate electrodes with positive f1 tuning and negative f2
            SID = ['EC' num2str(s)];
            idx = tbl.SID==s;
            data.resp = cat(1, data.resp, Dvow.(SID).resp(tbl.el(idx), :, :));
            respidx = [respidx; [repmat(s, sum(idx), 1), tbl.el(idx)]];
            clear idx SID
        end
        A = squeeze(mean(data.resp, 2));
        nancol = sum(isnan(A))>0.2*size(A, 1);
        A(:, nancol) = [];        
        nanrow = sum(isnan(A), 2)>0;

        % remaining electrodes after subsetting to electrodes with greatest
        % trial overlap
        remel = inflections(~nanrow, :);
        disp('-------------------------New elec set-----------------------------');
        disp(['Initial electrode set (' labels{ctr} '): ' num2str(size(data.resp, 1))]);    
        data.resp = data.resp(~nanrow, :, ~nancol);
        neur.respidx{ctr} = respidx(~nanrow, :);
        disp(['End electrode set (' labels{ctr} '): ' num2str(size(data.resp, 1))]);

        % run pairwise vowel decoding
        [svmOut{ctr}, confMat{ctr},~, ~, ~] = ecog_svm_main(data, 1:height(remel), ...
                                    'resp', 'vowel', 1:size(data.resp, 3), ...
                                    '', neur.tps, neur.pca, 1);

        subplot(3, 1, ctr)
        imagesc(max(confMat{ctr}, [], 3, 'omitnan')); 
        caxis([0.5 0.7]);  
        ctr = ctr+1;
        clear A nancol nanrow remel 
    end
    
    % populate the decoding structure    
    neur.confMat = confMat;
    neur.svmOut = svmOut;
    neur.bef = 20;
    neur.aft = 50;    
    save(filename, 'neur')
else
    load(filename)
end

types = {pf1nf2, nf1pf2, allsl, 1};

% revisualize accuracy as two clouds
agg = nan(10, 4);
y = nan(5, 4);
figure;
cols = getColors(1); 
for v = 1:5       
    for t = 1:3
        a = max(neur.confMat{t}, [], 3); 

        c = a;
        c(isnan(c))=0;
        c = c+c';
        c(c==0) = NaN;
        y(v, t) = mean(c(v, :), 'omitnan');
        
        a(isnan(a)) = [];
        agg(:, t) = a;
        if v==1
            text(t-0.05, mean(a)+0.025, ['n=' num2str(height(types{t}))]); 
            hold on;
        end
    end
    plot(1:4, y(v, :), 'Color', cols(v, :)); hold on;
    scatter(1:4, y(v, :), 85, cols(v, :), 'filled', 'MarkerFaceAlpha', 0.5); 
end   

% single electrode case
% for single electrode case
decode_model = ['decode_vowel_Spanish_' num2str(neur.tps(1)) '-' num2str(neur.tps(end))...
    neur.frametype '_allVowel_vowelftestEl-single_r0.15_.mat'];
a = load([datapath '/ecog_decode/' decode_model]);
neur_model{1} = a.neur;
acs_model{1} = a.acs;
sids = neur.respidx{1}(:, 1);
ssids = cellfun(@(x) str2double(x(3:end)), sSIDs);
idx = ismember(sids, ssids);
clear a
neur_model{1}.decode_single = table(sids(idx, 1), neur.respidx{1}(idx, 2), ...
    neur_model{1}.confMat(idx)', 'VariableNames',  {'SID', 'el', 'conf'});
          
decode_single = neur_model{1}.decode_single;
for i = 1:height(decode_single)
    conf = decode_single.conf{i};
    confSingle(:, :, i) = max(conf, [], 3, 'omitnan');   
end

confMax = nan(5, 5);
el = nan(5, 1);
for v = 1:5
    % maximum electrode per row
    confSingle(confSingle==0)= NaN;
    [~, el(v)] = max(mean(confSingle(v, :, :), 'omitnan'), [], 3, 'omitnan'); 
    confMax(v, :) = confSingle(v, :, el(v));
    
    c = confMax;
    c(isnan(c))=0;
    c = c+c';
    c(c==0) = NaN;
        
    plot([0 1], [mean(c(v, :), 'omitnan') y(v, 1)]', 'Color', cols(v, :)); hold on;
    scatter(0, (mean(c(v, :), 'omitnan')), 85, cols(v, :), ...
        'filled', 'MarkerFaceAlpha', 0.5); 
end

xticks(0:3);
xticklabels([{'single'}, labels]);
xlim([-0.5 3.5]);
ylim([0.5 0.66]);
set(gca, 'FontSize', 15);
yticks(0.5:0.05:0.7);

xlabel('Decoding Type');
ylabel('Pairwise Accuracy');
title('Decoding for different electrode types');
set(gca, 'FontSize', 15);
% plot(agg', 'o-', 'LineWidth', 0.5);


clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections* formant *model

%% Figure 4 - OPT: ERPs for single direction of encoding

neur = struct();
neur.tps = 20:50;
neur.frametype = 'all';
filename = [datapath ...
    '/ecog_decode/decode_vowel_Spanish_' num2str(neur.tps(1)) '-' num2str(neur.tps(end))...
    neur.frametype '_allVowel_slopeDirEl_r0.15_.mat'];
load(filename)

numvows = 5;
maxtp = nan(5, 1);
% find maximum time point per svm
svmOut = neur.svmOut{1};
classes = [svmOut.classesNum];
classes = mat2cell(classes, 2, ones(10, 1));
% comps = find(sum(classes==v-1)>0);

% length(tps)
allresp = nan(size(neur.respidx, 1), bef+aft+1, numvows, numvows);
for el = 1:size(decode_all.respidx, 1)
    confmat = nan(numvows, numvows, length(tps));
        
    % ADDED: multiply by high gamma amplitude
    sid = decode_all.respidx(el, 1);
    elnum = decode_all.respidx(el, 2);
    
    % translate electrode index from betas to decode single table
    sel =  find(decode_single.SID == sid & decode_single.el == elnum);
      
    % hg_weights = nan(numvows, length(tps));
    hg_weights = nan(numvows, bef+aft+1);
    % add resp to decode_single
    for v = 1:numvows
        SID = ['EC' num2str(sid)];
        Dvow = addtoDD(Dvow, 'dimex', bef, aft, {SID});
        hg_weights(v, :) = mean(Dvow.(SID).resp(elnum, ...
            :, Dvow.vowelType==v), 3, 'omitnan'); %tps     
    end
    
    for v = 1:numvows
        allresp(el, :, v, v) = hg_weights(v, :)'.* ...
            max(squeeze(decode_single.weights(sel, v, :))); % one time point
            %squeeze(decode_single.weights(sel, v, :));
        for v2=1:5
            if v ~= v2
                svmOut = decode_all.svmOut;
                
                % find the comparison with weights
                c = cellfun(@(x) all(x==[v-1; v2-1] | x==[v2-1; v-1]), classes);
                if all(classes{c} == [v-1; v2-1])
                    w = svmOut(c).beta(:, el); 
                else
                    w = -1*svmOut(c).beta(:, el); 
                end

                [~, idx] = max(abs(w));
                w = w(idx);
                allresp(el, :, v, v2) = hg_weights(v2, :)'.*w.*20;
                clear c w
            end
        end
    end
    
    decode_single.hg_resp(sel) = {hg_weights};        
end

% find maximum svm weight per electrode
Dvow.all.resp = nan(1, length(tps), numvows);

figure('Renderer', 'painters', 'Position', [200 800 200 350*4]);
% plot weighted average ERPs
cols = getColors(1);
for v = 1:numvows
    subplot(5, 1, v);
    
    for v2 = 1:numvows
        % x = (tps-bef)./100;
        x = (-bef:aft)./100;
        if v2==v
            y = squeeze(mean(allresp(:, :, v, v2), 'omitnan'));
            plot(x, y-y(1), 'LineWidth', 2.5, 'Color', cols(v2, :)); hold on;
        else
            plot(x, squeeze(mean(allresp(:, :, v, v2), 'omitnan')), 'LineWidth', 2, ...
                'Color', cols(v2, :), 'LineStyle', ':'); hold on;
        end
    end
    xticks([]);
    yticks([]);
    ylim([-0.1 0.45]);
    xlim([-0.2 0.5]);
    xline(0);
    yline(0);
    box off;
end
yticks(0:0.2:0.4);
xticks([-0.2 0 0.4]);
xlabel('Time (s)');
set(gca, 'FontSize', 15);
box off;

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections* formant *model

%% Figure 4 - OPT: Decoding from single electrodes vs. pooled across subject

% load in models

pitchSel = ''; % '', 'gt170', 'lt170'
decode_model = ['decode_vowel_Spanish_20-50_allVowel_vowelftestEl-single_r0.15_' ...
    pitchSel '.mat'];
load([datapath '/ecog_decode/' decode_model])

single = neur;
decode_model = ['decode_vowel_Spanish_20-50_allVowel_vowelftestEl-subj_r0.15_' ...
    pitchSel '.mat'];
load([datapath '/ecog_decode/' decode_model])
pooled = neur;
subjs = pooled.subj;

% for each subject
figure;
% 7 subjects
ctr = 1;
%subjs = {'EC100', 'EC214', 'EC172'};
for subj = subjs      
    subplot(2, 4, ctr);
    plot(squeeze(mean(pooled.confMat{ctr}, [1 2], 'omitnan')), ...
        'LineWidth', 2, 'Color', 'k');
    hold on;
    
    s = str2double(subj{1}(3:end));
    elecs = find(single.respidx(:, 1) == s);
    for e = elecs'
        plot(squeeze(mean(single.confMat{e}, [1 2],'omitnan')), ...
            'LineWidth', 0.75, ...
            'Color', [0.5 0.5 0.5]); hold on;
    end
    ylim([0.4 0.7]);
    yline(0.5, 'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'r');
    title(subj{1});
    
    ctr = ctr + 1;
end
sgtitle(pitchSel)

% new visualization in confusion matrix format
figure;
% 7 subjects
comp = [2 3 4 5 8 9 10 14 15 20];
pairs = [1, 2; 1, 3; 1, 4; 1, 5; 2, 3; 2, 4; 2, 5; 3, 4; 3, 5; 4, 5];
for subj = subjs  
    s = str2double(subj{1}(3:end));
    elecs = find(single.respidx(:, 1) == s);
     
    % find max value per comparison for elecs in this subject
    tmp = cell2mat(cellfun(@(x) max(x, [], 3, 'omitnan'), single.confMat(elecs), ...
        'UniformOutput', false));
    sconf = max(reshape(tmp, [5, 5, length(elecs)]), [], 3, 'omitnan');
    
    % get max tp confMat for pooled electrodes in this subject
    sidx = find(strcmp(pooled.subj, subj));
    pconf = max(pooled.confMat{sidx}, [], 3, 'omitnan');
    
    for ctr = 1:10      
        subplot(4, 5, comp(ctr));
                
        x1 = 1-0.1+rand*0.2;
        x2 = 2-0.1+rand*0.2;
        y1 = sconf(pairs(ctr, 1), pairs(ctr, 2));
        y2 = pconf(pairs(ctr, 1), pairs(ctr, 2));
        scatter(x1, y1, 35, ...
            [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.5); hold on;
        scatter(x2, y2, 35, ...
            'k', 'filled', 'MarkerFaceAlpha', 0.5); hold on;
        plot([x1 x2], [y1 y2], 'LineWidth', 0.2, 'Color', 'k');
        
        ylim([0.45 0.8]);
        xlim([0.5 2.5]);
        yline(0.5, 'LineWidth', 1.5, 'LineStyle', '--', 'Color', 'r');
    end
    sgtitle(pitchSel)
end

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections* formant *model

%% Figure 4 - OPT: Weighted Average of Comparison HG RF using SVM weights

load('out_elecs_voweltypeftest_bychan.mat')

decode_single = neur_model{1}.decode_single;
decode_all = neur_model{2};
tps = neur_model{1, 1}.tps;

% each svm use only one time point
weights = nan(height(decode_single), 10);
keep = zeros(height(decode_single), 10);
hg = zeros(height(decode_single), 10);
imgs = cell(height(decode_single), 10);

numvows = 5;

[N, XEDGES, YEDGES, BINX, BINY] = getF12Bin(Dvow);

% get weights per electrode
svmOut = decode_all.svmOut;
for c = 1:size(svmOut, 2)
        
    % find maximum timepoint for this SVM
    [~, maxidx] = max(mean(svmOut(2).acc, 2));
    
    % find the weight per electrode
    for el = 1:size(decode_all.respidx, 1)
        confmat = nan(numvows, numvows, length(tps));
        
        w = svmOut(c).beta(maxidx, el);
        v1 = svmOut(c).classesNum(1) + 1;
        v2 = svmOut(c).classesNum(2) + 1;
        confmat(v1, v2, :) = w;
        confmat(v1, v2, :) = -w;             

        % translate electrode index from betas to decode single table
        sid = decode_all.respidx(el, 1);
        elnum = decode_all.respidx(el, 2);      
        sel =  find(decode_single.SID == sid & decode_single.el == elnum);
        weights(sel, c) = w*30;
        
        % decide whether to keep this electrode
        keep(sel, c) = w>prctile(svmOut(c).beta(maxidx, :), 90) | ...
             w<prctile(svmOut(c).beta(maxidx, :), 10);
        
        
        % get all electrode RF
        if keep(sel, c)
            SID = ['EC' num2str(sid)];
            bef = 20;
            hgmean = squeeze(mean(Dvow.(SID).resp(elnum, ...
                maxidx+bef-2:maxidx+bef+2, :), 2, 'omitnan'));
            img=nan(size(N, 1),size(N, 2));
            for y=1:length(YEDGES)-1
                for x=1:length(XEDGES)-1
                    if sum(BINX==x & BINY==y, 'omitnan')>5               
                        % mean hg peak as pixel value
                        img(x, y)=mean(hgmean(BINX==x & BINY==y), 'omitnan');
                    end
                end
            end
            hg(sel, c) = mean(hgmean, 'omitnan');
            imgs{sel, c} = smoothdata(img);           
            
        end
    end
end
sum(keep)

% weight RFs and average across electrodes, plot
pos = reshape(1:25, 5, 5)';
figure;
medoids = plotVariance(Dvow.formantVals, Dvow.vowel, getColors(1), 0, 0); 
ys=arrayfun(@(x) find(abs(x-YEDGES)==min(abs(x-YEDGES))), medoids(:, 1));
xs=arrayfun(@(x) find(abs(x-XEDGES)==min(abs(x-XEDGES))), medoids(:, 2));
bm = getColors(1);
avgs = cell(5, 5);
for col = 1:10
    
    v1 = svmOut(col).classesNum(1) + 1;
    v2 = svmOut(col).classesNum(2) + 1;
        
    subplot(4, 5, pos(v1, v2));
    avg=zeros(size(N, 2),size(N, 1));
    ridx = ones(height(imgs), 1);
    for row = find(keep(:, col))'
        t = imgs(row, col);
        m = (t{1}*weights(row, col))';
        avg = avg + m*hg(row, col)*15;
        avg(isnan(avg)) = m(isnan(avg));
    end
    
    imagesc(avg); hold on;
    avgs{v1, v2}=avg;
    
    maxacc = max(mean(svmOut(col).acc, 2));
    title([svmOut(col).classes{1} '-' svmOut(col).classes{2} ': ' num2str(maxacc)]);
    set(gca,'XDir','reverse');
    colormap( [1 1 1; flipud(pink(256))]);  
    
    scatter(xs, ys, 25, bm, 'filled')
    scatter(xs(v1), ys(v1), 105, bm(v1, :), 'filled', 'MarkerEdgeColor', 'k');
    scatter(xs(v2), ys(v2), 105, bm(v2, :), 'filled', 'MarkerEdgeColor', 'k');
end

% average for each vowel
fig1=figure('Position', [100 800 300*5 250]);
fig2=figure('Position', [100 800 300*5 250]);
vowavg=cell(5, 1);
% a e i o u
% a-i, a-u (1-3, 1-5)
exclude = [];
%avgs{1, 5} = [];
avgs{1, 3} = [];
for v = 1:5
   pos = avgs(:, v);
   neg = avgs(v, :);
   
   for y = 1:5
       % if adjacent
       if ~isempty(pos{y})
           if isempty(vowavg{v})
               vowavg{v} = pos{y};
           else
               vowavg{v} = vowavg{v} + pos{y};
           end
       end
       
       if ~isempty(neg{y})
           if isempty(vowavg{v})
               vowavg{v} = neg{y};
           else
               vowavg{v} = vowavg{v} - neg{y};
           end          
       end       
   end
   set(0, 'CurrentFigure', fig1);
   subplot(1, 5, v);
   avg = vowavg{v};
   imagesc(XEDGES./1000, YEDGES./1000, smoothdata(vowavg{v})); hold on;      
   set(gca,'XDir','reverse');
   colormap( [1 1 1; flipud(pink(256))]);  

   scatter(medoids(:, 2)./1000, medoids(:, 1)./1000, 25, bm, 'filled');
   scatter(medoids(v, 2)./1000, medoids(v, 1)./1000, 105, bm(v, :), 'filled', 'MarkerEdgeColor', 'k');
   box off;
   xticks([]);
   yticks([]);
   
   % plot the convex hull
   frm = [Dvow.formantVals(1, Dvow.vowelType==v); Dvow.formantVals(2, Dvow.vowelType==v)];
   frm = rmoutliers(frm', 'quartiles', 'ThresholdFactor', 0.85); %'percentiles',[5 95]
   [k,av] = convhull(frm);
   plot(frm(k, 2)./1000, frm(k, 1)./1000, 'Color', bm(v, :), 'LineWidth', 2);
   
   set(0, 'CurrentFigure', fig2);
   subplot(2, 5, v+5);
   d1 = sum(avg, 'omitnan') ./ sum(~isnan(avg), 'omitnan');
   plot(XEDGES(1:end-1)+diff(XEDGES)/2, d1, 'LineWidth', 1.5, 'Color', 'k');
   d1_ac = histcounts(discretize(Dvow.formantVals(2,...
        Dvow.vowelType == v), XEDGES), 0:length(XEDGES));
   yyaxis right
   plot(XEDGES, d1_ac, 'LineWidth', 1.5, 'Color', bm(v, :));
   xlim([600 3000])

   d2 = (sum(avg, 2, 'omitnan') ./ sum(~isnan(avg), 2, 'omitnan'))';
   subplot(2, 5, v);
   plot(YEDGES(1:end-1)+diff(YEDGES)/2, d2, 'LineWidth', 1.5, 'Color', 'k');
   yyaxis right
   d2_ac = histcounts(discretize(Dvow.formantVals(1,...
       Dvow.vowelType == v), YEDGES), 0:length(YEDGES));
   plot(YEDGES, d2_ac, 'LineWidth', 1.5, 'Color', bm(v, :));
   xlim([200 850])
   clear d*       
end

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections* formant *model

%% Figure 4 - OPT: Weighted Average of Vowel HG RF using SVM weights

load('out_elecs_voweltypeftest_bychan.mat')

decode_single = neur_model{1}.decode_single;
decode_all = neur_model{2};
tps = neur_model{1, 1}.tps;
weights = nan(height(decode_single), 5, length(tps));
debug = 0;
numvows = 5;
for el = 1:size(decode_all.respidx, 1)
    svmOut = decode_all.svmOut;
    confmat = nan(numvows, numvows, length(tps));
    for c = 1:size(svmOut, 2)
        w = svmOut(c).beta(:, el); 
        x = svmOut(c).classesNum(1) + 1;
        y = svmOut(c).classesNum(2) + 1;
        confmat(x, y, :) = w;
        confmat(y, x, :) = -w;
    end
    
    % ADDED: multiply by high gamma amplitude
    sid = decode_all.respidx(el, 1);
    
    elnum = decode_all.respidx(el, 2);
    hg_weights = nan(numvows, length(tps));
    for v = 1:numvows
        SID = ['EC' num2str(sid)];
        Dvow = addtoDD(Dvow, 'dimex', bef, aft, {SID});
        hg_weights(v, :) = mean(Dvow.(SID).resp(elnum, ...
            tps, Dvow.vowelType==v), 3, 'omitnan');
    end
    
    % translate electrode index from betas to decode single table
    sel =  find(decode_single.SID == sid & decode_single.el == elnum);
    weights(sel, :, :) = squeeze(sum(confmat, 'omitnan')).*(20*hg_weights);
    
    if debug
        f = figure;
        subplot(1, 3, 1);
        imagesc(corr(squeeze(weights(sel, :, :)),hg_weights));
        
        subplot(1, 3, 2);
        cols = getColors(1);
        for v = 1:5            
            plot(squeeze(weights(sel, v, :))', 'Color', cols(v, :)); hold on;
        end
        legend({'a', 'e', 'i', 'o', 'u'});
        
        subplot(1, 3, 3);
        plotVowelErp(Dvow, ['EC' num2str(decode_single{sel, 'SID'})], ...
             decode_single{sel, 'el'}, 'dimex', f)
    end
end

decode_single.weights = weights;
decode_single.tps = repmat(tps, height(decode_single), 1);
neur_model{1, 1}.decode_single = decode_single;

mSIDs = arrayfun(@(x) ['EC' num2str(x)], ...
    unique(decode_all.respidx(:, 1)), 'UniformOutput', false);
imagescVowelResponse(Dvow, mSIDs', allidx, unique(Dvow.vowelType), decode_single);
clear mSIDs

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections* formant *model



%% ------------------------ Functions -------------------------------------

function [allidx, fvals] = saveElecs(SIDs, corpus, Dvow)
    % finds subject specific max time point
    info=struct();
    info.resp='resp'; 
    [allidx, fvals] = getElecs(Dvow, SIDs, [], corpus, 'ftest', info);
    save('out_elecs_voweltypeftest_bychan.mat', 'allidx', 'fvals');
end

function neur = lda(A, pcaflag, idx, Svow)

    max_comp = 800;

    % remove columns and rows where data is missing
    % trials for which over electrodes do not have data
    % remanining electrodes with missing trials
    nancol = sum(isnan(A))>0.35*size(A, 1);
    A(:, nancol) = [];
    nanrow = sum(isnan(A), 2)>0;
    A(nanrow, :) = [];
    
    % find which trials retained
    tmp = find(idx);
    trls = tmp(~nancol);
    
    % populate neural LDA model structure
    neur = struct();
    neur.trls = trls;
    neur.y = Svow.vowelType(trls)';
    
    if pcaflag
        % with pca
        expl_perc = 90;
        [neur.coeff, ~, ~, ~, exp] = pca(A, 'NumComponents', max_comp); 
        neur.n_comp = find(diff(cumsum(exp)>expl_perc)); %size(neur.coeff, 2); %
        clear score exp
    else
        neur.coeff = A;
        neur.n_comp = size(A, 2);
    end
    
    % run lda with kfold cross-validation
    neur.kfolds = 5;
    [neur.y_pred, neur.mappedX, ~] = LDAmap(neur.coeff(:, 1:neur.n_comp), ...
        neur.y, neur.kfolds);
    neur.conf = squareform(pdist(squeeze(neur.mappedX(1, :, 1:2))));
end

% find the subplot number based on a grid formant of the subplot
function plotNum = gridCount(row, col, ~, cols)
    plotNum =  (row-1)*cols + col;
end
