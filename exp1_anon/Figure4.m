
% "STARTUP.M" NEEDS TO BE RUN BEFORE THE FOLLOWING CODE

%% ------------ Vowel Formant Encoding in Natural Speech ------------------

[inflections]=loadInfl(betaInfo, SIDs, 1, Dvow, bef, aft);

% '', 'gt170', 'lt170'
pitchSel = ['' '_anon']; 

decode_model = ['decode_vowel_Spanish_20-50_allVowel_vowelftestEl-single_r0.15_' pitchSel '.mat'];
load([datapath '/ecog_decode/' decode_model])

neur_model{1} = neur;
acs_model{1} = acs;
          
sids = neur.respidx(:, 1);
els = neur.respidx(:, 2);

ssids = cellfun(@(x) str2double(x(2:end)), SIDs);
idx = ismember(sids, ssids);

neur_model{1}.decode_single = table(sids(idx), els(idx), neur.confMat(idx)', ...
    'VariableNames',  {'SID', 'el', 'conf'});

% for pooled within subj: 'decode_vowel_Spanish_20-50_allVowel_vowelftestEl-subj_r0.15_'
% only 5 subjects used to maximize overlapping trials
% for pooled across subjs: 'decode_vowel_Spanish-mini_20-50_allVowel_vowelftestEl_r0.15_'
decode_model = ['decode_vowel_Spanish-mini2_20-50_allVowel_vowelftestEl_r0.15_' pitchSel '.mat'];
load([datapath '/ecog_decode/' decode_model])

% clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* inflections ...
%     formant *model 

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

% single electrode decoding
filename_single = ...
    'ecog_decode/decode_vowel_Spanish_maxFstattp_allVowel_vowelftestEl-single_.mat';
% population of electrodes decoding
filename_pop = ...
    'ecog_decode/decode_vowel_Spanish_maxFstattp-pca_allVowel_vowelftestEl-both.mat';

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
    boxchart(ones(50, 1)*v, squeeze(pairavg(v, :))', 'MarkerColor', ...
        'none', 'LineWidth', 0.5, ...
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
disp(['Pairs at p<.0001 = ' num2str(sum(p<0.0001)) '/10']);
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

    SID = ['S' num2str(inflections.SID(els(ctr)))];
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
elidx.SID = arrayfun(@(x) ['S' num2str(x)], elidx.SID, ...
    'UniformOutput', false);

vals = plotAvgImagescResponse(Dvow, elidx, getColors(1), [], ...
    rescale(max([elidx.sr_f1, elidx.sr_f2], [], 2)));
title('F1-/F2+');
set(gca, 'FontSize', 15);
caxis([-0.0 0.45])
% cmap = brewermap(256, 'PuOr'); 
% colormap([1 1 1; flipud(cmap(40:240, :))]);
box off;
yticks([]);
xticks([]);
xlim([800, 3000]);
ylim([200 800]);
brighten(0.1);

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
elidx.SID = arrayfun(@(x) ['S' num2str(x)], elidx.SID, ...
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
    [h,p, ~, stats] = ttest2(Y(1, :), Y(2, :));
    pair_label = strjoin(decode_all.class_ord(i).classes);
    supp_stats = [supp_stats; {pair_label, mean(Y, 2), stats.tstat, p}];
end
disp(supp_stats)

load('out_elecs_voweltypeftest_bychan_anon.mat')

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

% find the subplot number based on a grid formant of the subplot
function plotNum = gridCount(row, col, ~, cols)
    plotNum =  (row-1)*cols + col;
end
