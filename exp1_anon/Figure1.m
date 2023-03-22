
% "STARTUP.M" NEEDS TO BE RUN BEFORE THE FOLLOWING CODE

%% Figure 1 (a, b) - Single Sentence Vowel Formants

% setup
cols = getColors(1);
vows = {'a', 'e', 'i', 'o', 'u'};

load([datapath 'stim_info/mel_centerF.mat']);
aud = binfrqs(1:end-1)+diff(binfrqs)/2;
aud = aud(1:end-1);

out = load([datapath 'stim_info/out_sentence_details_dimex_all_loudness.mat']);
sent = 25;
x = out.sentdet(sent);

width = 455*2;
height = 200*2;
figure('Position', [100 800 width height]);
text(0,0.5,{'fondo de población de las naciones unidas'}, 'FontSize', 15);
set(gca,'Visible','off')

% sound wave form
subplot(3, 3, [1 2]);
soundf = x.soundf;
ons = x.soundOns;
befaft = x.befaft;
range = round([(befaft(1)+ons(1))*soundf ...
    (befaft(1)+ons(2))*soundf]);
plot(range(1)/soundf:1/soundf:range(2)/soundf, x.sound(range(1):range(2)), ...
    'Color', [0.7, 0.7, 0.7]);

transf = soundf / x.dataf;
vowtimes = x.vowelTimes;
for i=1:size(x.vowelTimes, 2)
    col=find(contains(vows, x.vowel{i}));
    %xval = ceil(vowtimes(1, i)*transf)-range(1);
    xval = vowtimes(1, i)/x.dataf;
    xline(xval, 'Color', cols(col, :), 'LineWidth', 0.5);
    text(xval-500, 0.6, x.vowel{i}, ...
        'FontWeight', 'bold', 'FontSize', 15, 'Color', cols(col, :));
end
set(gca,'Visible','off')
xlim([find(x.onsOff(1, :))/100 ...
    find(x.onsOff(2, :))/100] );

% spectrogram of the sentence
subplot(3, 3, [4 5]);
imagesc((1:size(x.aud, 2))/100, aud./1000, x.aud);
colormap(flipud(gray))
for i=1:size(x.vowelTimes, 2)
    col=contains(vows, x.vowel{i});
    xline(x.vowelTimes(1, i)/100, ...
        'Color', cols(col, :), 'LineWidth', 2);
end
xlim([find(x.onsOff(1, :))/100 ...
    find(x.onsOff(2, :))/100] );
set(gca, 'YDir', 'normal')
xticks([]);
ylim([0 3]);
yticks(0:1:3)

freq_ranges = [prctile(Dvow.formantVals(1:2, :)', 2.5) ...
    prctile(Dvow.formantVals(1:2, :)', 97.5)]';
horzline(freq_ranges./1000, gcf, 'r', '-', 1);

ylabel({'Frequency', '(kHz)'});

box off
set(gca, 'FontSize', 13)

% Figure 0 - Formants for Single Sentence
subplot(3, 3, [7 8]);
flim = [700 2200];
fcols =  [0.1, 0.1, 0.1; 0.5, 0.5, 0.5];
for f=1:2
    for i=1:size(x.vowelTimes, 2)
        col=contains(vows, x.vowel{i});

        vowelRange = x.vowelTimes(1, i):x.vowelTimes(2, i);
        plot(vowelRange/100, x.formants(f, vowelRange), ...
            'Color', cols(col, :), 'LineWidth', 1.5); hold on;
        
        medtp = ceil(median(1:length(vowelRange)));
        scatter(vowelRange(medtp)/100, x.formants(f, vowelRange(medtp)), ...
            35, fcols(f, :), 'filled');
    end
end
ylabel({'Frequency', '(kHz)'});
box off
yline(700)
yticks([1000 2000]); yticklabels({'1', '2'})
set(gca, 'FontSize', 13)
xlim([find(x.onsOff(1, :))/100 ...
    find(x.onsOff(2, :))/100] )
xlabel('Time (s)');

% Figure 0 - Spanish Corpus Formant Variation
subplot(3, 3, [3 6 9])
plotVariance(Dvow.formantVals, Dvow.vowel, cols, 0, 1, 1);
set(gca, 'FontSize', 15);
xticks(1000:1000:3000);
xticklabels(split(num2str(1:1:3)));
yticks(200:400:1000);
yticklabels(split(num2str(0.2:0.4:1)));
ylabel('F1 (kHz)');
xlabel('F2 (kHz)');
legend('Location', 'southeast');

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* inflections ...
    formant *model

%% Figure 1 (c) - Spanish vowel spectrum

cols = getColors(1);
vows = {'a', 'e', 'i', 'o', 'u'};
load([datapath 'stim_info/mel_centerF.mat']);

fs = Dvow.meanf0>170;
x = binfrqs(1:end-1)+diff(binfrqs)/2;
x = x(1:end-1);
figure('Position', [800 200 100 700]);
for v = 1:5
    subplot(5, 1, v);
    title(vows{v})
    ms = Dvow.meanf0<170;

    % remove outliers on the sample dimension
    cl = rmoutliers(squeeze(median(Dvow.aud(:, 50:55, ...
        strcmp(Dvow.vowel, vows{v}) & ms), 2))');

    % pre-emphasis accounting for 1/f fall-off 
    % (based on https://www.fon.hum.uva.nl/praat/manual/Sound__Filter__pre-emphasis____.html)
    alpha  = exp(-2*pi*3*(1/100));

    cl = cl(:, 1:79)-alpha*cl(:, 2:80);

    plot(x(2:end)./1000, smoothdata(median(cl), 'SmoothingFactor',0.2), ...
        'Color', cols(v, :), 'LineWidth', 2);

    % showing location of formant median
    title(vows{v})
    xlim([0.25, 2.5]);
    box off;
    if v<5, xticks([]); end
end

[dim_r, dim_p] = corr(Dvow.formantVals(1, :)', Dvow.formantVals(2, :)');
[tim_r, tim_p] = corr(Tvow.formantVals(1, :)', Tvow.formantVals(2, :)');
disp(['Dimex F1-F2 correlation:' num2str(dim_r) ' and p-val: ' num2str(dim_p)]);
disp(['Timit F1-F2 correlation:' num2str(tim_r) ' and p-val: ' num2str(tim_p)]);

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* inflections ...
    formant *model

%% Figure 1 (d, e, f) - Example Vowel ERPs
addpath(genpath('util'))
Dvow = addtoDD(Dvow, 'dimex', bef, aft, {'S1', 'S7'});
cols = getColors(1);

% E1
erpSIDs = {'S1', 'S1', 'S1'}; 
els = [55 71, 118];

fig = figure('Position', [100, 600, 800, 200]);
for i = 1:3
    el = els(i);
    SID = erpSIDs{i};

    subplot(1, 3, i)
    Dvow.(SID).resp(el, :, :) = normalize(Dvow.(SID).resp(el, :, :));
    
    plotVowelErp(Dvow, SID, el, fig, cols, bef/100);
    ylabel('HGA (z)')
    set(gca, 'FontSize', 15);
    yticks(-0.5:0.5:1);
    title(['E' num2str(i)]);
    legend('off')
end

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* inflections ...
    formant *model desel


%% Figure 1 (g) - Vowel F-Stat on Native Brain
elecs = containers.Map;
elecs('S1') = [55 71 118]; % E1 and E2, ERP electrode  118
elecs('S7') = 54; % E3, Sigmoid 54

% initialize design electrode structure
desel=struct();
numbins = 30;
desel.conds = 1:numbins;

load([datapath 'pt_data/out_elecs_speechtypeftest_bychan_dimex_Spanish_anon.mat'], ...
    'allidx');
allidx_speech = allidx;

load('pt_data/out_elecs_voweltypeftest_bychan_anon.mat', 'allidx', 'fvals');
allidx_vowel = allidx;

modelname = {'onset_maxDtL_maxDtLOnset_vowelOnset_aud'};
for s = SIDs
    SID = s{1};
    strfs = loadMultModelStrf(SID, modelname, 'dimex', datapath, 1, '');
end 

[stgelecs,~] = getElecs(Dvow, SIDs, 'dimex', 'anatomy', struct(), ...
    [datapath '/pt_data']);

varnames = {'SID', 'el', 'speech', 'stg'};
subelecs = array2table(zeros(0,4), 'VariableNames', varnames);
speechelecs = array2table(zeros(0,4), 'VariableNames', varnames);
for s = SIDs
    SID = s{1};
    sids = repmat(str2double(SID(2:end)), length(allidx_vowel.(SID)), 1);
    t2 = table(sids, allidx_vowel.(SID)', ismember(allidx_vowel.(SID)', ...
        allidx_speech.(SID)), ismember(allidx_vowel.(SID), stgelecs.(SID))', ...
        'VariableNames', varnames);
    subelecs = [subelecs; t2];
    
    sids = repmat(str2double(SID(2:end)), length(allidx_speech.(SID)), 1);
    t2 = table(sids, allidx_speech.(SID), ismember(allidx_speech.(SID), ...
        allidx.(SID)), ismember(allidx_speech.(SID), stgelecs.(SID)), ...
        'VariableNames', varnames);
    speechelecs = [speechelecs; t2];
    clear t2
end
clear t2 s

disp('----------------- Subject Stats -------------------------')
disp(['Total vowel electrodes: ' num2str(height(subelecs))]);
disp(['Total subj: ' num2str(length(unique(subelecs.SID)))]);
disp(['Min vowel electrodes per subj: ' num2str(min(crosstab(subelecs.SID)))]);
disp(['Max vowel electrodes per subj: ' num2str(max(crosstab(subelecs.SID)))]);
disp(['Median vowel electrodes per subj: ' num2str(median(crosstab(subelecs.SID)))]);
disp('------------- Speech Stats -------------')
disp(['Total speech responsive: ' num2str(height(speechelecs))]);
disp(['Total speech responsive on STG: ' num2str(sum(speechelecs.stg))]);
disp(['Count vowel responsive: ' num2str(height(subelecs))]);
disp(['Count vowel responsive on STG: ' num2str(sum(subelecs.stg))]);

tmp = linspace(6.5, 30, numbins-2);
edges = [0 tmp 400];
for s=1:length(SIDs)
    SID = SIDs{s};
    if isfield(betaInfo, SID) && ~isempty(betaInfo.(SID).els)
        conds = discretize(fvals.(SID), edges); 
        conds(1, conds>1) = conds(1, conds>1)+1;
        
        % the second condition corresponds to speech responsive electrodes
        conds(1, intersect(find(conds==1),allidx_speech.(SID))) = 2;
        desel.(SID).elid = 1:length(fvals.(SID));
        desel.(SID).condition=conds;
    end        
end

for k = elecs.keys()
    % selected electrodes
    key = k{1};
    desel.(key).selid = elecs(key);
end

% for single subject
desel.sz = [5, 5, repmat(55, 1, numbins)];
desel.cols = [ 0 0 0; 0 0 0 ; brewermap(numbins, 'YlGn')];
[~] = plotNativeElec({'S1'}, desel, [datapath '/pt_data']);
% print(fullfile('S1_vowel.jpg'), '-djpeg', '-painters', '-r600')
set(gcf,'Color','w');
set(gcf, 'InvertHardcopy', 'off');

% for the colorbar
figure;
colormap(desel.cols(2:end, :))
h=colorbar;
set(gca,'ColorScale','log');
caxis([edges(2) edges(end-1)]);
ylabel(h, 'Vowel F-stat');
set(gca, 'FontSize', 15);
set(h,'XTick',[5 30]);


clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections formant *model desel

%% Figure 1 - SUPP: Formant TRF model to vowel responsive
% run section above before running this section 

elecs = containers.Map;
elecs('S1') = []; % E1 - E3, ERP electrodes 55 71 118
% elecs('S7') = 54; % E3, Sigmoid

typeinfo = struct();
typeinfo.resp = 'resp';

typeinfo.modelname = {'onset_maxDtL_maxDtLOnset_vowelOnset_F0_audF1'};
[~, rsqs_f1] = getElecs(Dvow, SIDs, 'dimex', 'strfrsq', typeinfo, ...
    [datapath '/pt_data']);
typeinfo.modelname = {'onset_maxDtL_maxDtLOnset_vowelOnset_F0_audF2'};
[~, rsqs_f2] = getElecs(Dvow, SIDs, 'dimex', 'strfrsq', typeinfo, ...
    [datapath '/pt_data']);
typeinfo.modelname = {'onset_maxDtL_maxDtLOnset_vowelOnset_F0_audF1F2'};
[~, rsqs] = getElecs(Dvow, SIDs, 'dimex', 'strfrsq', typeinfo, ...
    [datapath '/pt_data']);
typeinfo.modelname = {'onset_maxDtL_maxDtLOnset_vowelOnset'}; %_vowelOnset_F0
[~, rsqs_subform] = getElecs(Dvow, SIDs, 'dimex', 'strfrsq', typeinfo, ...
    [datapath '/pt_data']);
[stgelecs,~] = getElecs(Dvow, SIDs, 'dimex', 'anatomy', struct(), ...
    [datapath '/pt_data']);
load('out_elecs_voweltypeftest_bychan_anon.mat', 'fvals', 'allidx');

figure('Position', [100 800 300 300], 'Renderer', 'Painters'); 
fvals_comp=[];
rsq_comp =[];
for s = SIDs
    SID = s{1};
    if isfield(rsqs_subform, SID)
        %cols = [0.5 0.8 0.6; 0.8 0.15 0.95];
        stg = ismember(1:length(fvals.(SID)), stgelecs.(SID));
        % take the maximum unique variance from each of the three models (F1, F2, F1+F2)
        % to account for single and co-encoding electrodes
        rsq = rsqs.(SID);
        
        selidx = allidx.(SID);
        idx = stg & fvals.(SID)>0.2 & ~ismember(1:length(fvals.(SID)), selidx);      
                
        scatter3(fvals.(SID)(idx), rsq(idx), find(idx), 45, [0.5 0.5 0.5], 'filled', ...
            'MarkerFaceAlpha', 0.7, 'LineWidth', 0.025); hold on;
        %text(fvals.(SID)(idx), rsq(idx), repmat({SID}, 1, sum(idx)))

        scatter3(fvals.(SID)(selidx), rsq(selidx), selidx, 45, [0.2 0.65 0.3], 'filled', ...
            'MarkerFaceAlpha', 0.7, 'LineWidth', 0.025); hold on;
        view(2)
        fvals_comp = [fvals_comp fvals.(SID)(selidx)];
        rsq_comp = [rsq_comp rsq(selidx)];
        
        if ismember(SID, elecs.keys)
            for e = elecs(SID)
                scatter3(fvals.(SID)(e), rsq(e), e, 55, [0.9 0 0], ...
                    'LineWidth', 1.2, 'MarkerEdgeColor', 'r'); hold on;
                text(fvals.(SID)(e), rsq(e), [SID ' e:' num2str(e)]);
            end
        end
    end
end

ylabel('Spectral TRF R^2');
xlabel('Vowel F-Stat');
grid off;
set(gca, 'Xscale', 'log', 'FontSize', 15);
[r, pval]= corr(fvals_comp',rsq_comp','Type', 'Spearman');
disp(['Rho: ' num2str(r) ', pval: ' num2str(pval)]);

desel.sz = [NaN, 3, repmat(55, 1, length(desel.conds))];
desel.cols = [ 0 0 0; 0 0 0; brewermap(length(desel.conds), 'YlGn')];
desel.labels = split(num2str(1:length(desel.conds)));
desel.EC100.selid = [];
desel.EC214.selid = [];

% colored by distance from the diagonal
pstg_thresh = -6;
[lh_mni] = plotMNIElec(SIDs, desel, 'lh', [datapath '/pt_data']);

% splitting pSTG and a-mSTG
% idx = lh_mni.cond>2 & lh_mni.y>=pstg_thresh;
% scatter3(lh_mni.x(idx), lh_mni.y(idx), lh_mni.z(idx), 35, 'r');

% plot3([-250, -265], [-30 -48], [4 2], 'LineWidth', 2, 'Color', 'r')
axes('Position',[.7 .1 .3 .3])
p = pie([sum(lh_mni.cond>2), sum(lh_mni.cond==2)], [1 1]); 
p(1).FaceColor = desel.cols(20, :);
p(1).EdgeColor = 'none';
p(3).FaceColor = [0.6 0.6 0.6];
p(3).EdgeColor = 'none';
% print(fullfile('mniBrain_LH.jpg'), '-djpeg', '-painters', '-r600');

[rh_mni] = plotMNIElec(SIDs, desel, 'rh', [datapath '/pt_data']);

% splitting pSTG and a-mSTG
% idx = rh_mni.cond>2 & rh_mni.y>=pstg_thresh;
% scatter3(rh_mni.x(idx), rh_mni.y(idx), rh_mni.z(idx), 35, 'r');

axes('Position',[.1 .1 .3 .3])
p = pie([sum(rh_mni.cond>2), sum(rh_mni.cond==2)], [1 1]); 
p(1).FaceColor = desel.cols(20, :);
p(1).EdgeColor = 'none';
p(3).FaceColor = [0.6 0.6 0.6];
p(3).EdgeColor = 'none';
% print(fullfile('mniBrain_RH.jpg'), '-djpeg', '-painters', '-r600')

figure;
mnis = {lh_mni, rh_mni};
titles = {'Left Hemi', 'Right Hemi'};
for hemi = 1:2
    subplot(1, 2, hemi)
    y = mnis{hemi}.y;

    %ys = rescale([y(mnis{hemi}.cond>2); pstg_thresh]);
    ys{hemi} = [y(mnis{hemi}.cond>2 & ...
        strcmp(mnis{hemi}.anatomy, 'superiortemporal'))];
    histogram(ys{hemi}, 6, 'FaceAlpha', 0.5, ...
        'FaceColor', desel.cols(20, :), 'EdgeColor', 'none', ...
        'Normalization', 'Probability'); % 

    xline(pstg_thresh, 'LineWidth', 2);
    xticks([-20 20]);
    xlim([-20 20]);

    ylabel('Probability')
    box off;
    xticklabels({'posterior', 'anterior'});
    set(gca, 'FontSize', 15, 'Xdir', 'reverse');
    title(titles{hemi});
end


disp('--------------- Counts -------------------')
disp(['LH vowel electrode count: ' num2str(sum(lh_mni.cond>2)) ...
    ' RH count: ' num2str(sum(rh_mni.cond>2))]);

disp('--------------- Stats -------------------')
[h, p] = kstest2(ys{1}, ys{2});
disp(['P-value for LH versus RH y-coordinate distribution: ' num2str(p)]);


clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections formant *model desel

%% Figure 1 - SUPP: All subject speech responsive vs. vowel responsive

% initialize design electrode structure
desel=struct();
numbins = 30;
desel.conds = 1:numbins;

load([datapath 'pt_data/out_elecs_speechtypeftest_bychan_dimex_Spanish_anon.mat'], ...
    'allidx', 'fvals');
allidx_speech = allidx;
fvals_speech = fvals;

load('pt_data/out_elecs_voweltypeftest_bychan_anon.mat', 'allidx', 'fvals');
allidx_vowel = allidx;
fvals_vowel = fvals;

modelname = {'onset_maxDtL_maxDtLOnset_vowelOnset_aud'};
for s = SIDs
    SID = s{1};
    strfs = loadMultModelStrf(SID, modelname, 'dimex', datapath, 1, '');
end 

[stgelecs,~] = getElecs(Dvow, SIDs, 'dimex', 'anatomy', struct(), ...
    [datapath '/pt_data']);

varnames = {'SID', 'el', 'speech', 'stg'};
subelecs = array2table(zeros(0,4), 'VariableNames', varnames);
speechelecs = array2table(zeros(0,4), 'VariableNames', varnames);
for s = SIDs
    SID = s{1};
    sids = repmat(str2double(SID(3:end)), length(allidx_vowel.(SID)), 1);
    t2 = table(sids, allidx_vowel.(SID)', ismember(allidx_vowel.(SID)', ...
        allidx_speech.(SID)), ismember(allidx_vowel.(SID), stgelecs.(SID))', ...
        'VariableNames', varnames);
    subelecs = [subelecs; t2];
    
    sids = repmat(str2double(SID(3:end)), length(allidx_speech.(SID)), 1);
    t2 = table(sids, allidx_speech.(SID), ismember(allidx_speech.(SID), ...
        allidx.(SID)), ismember(allidx_speech.(SID), stgelecs.(SID)), ...
        'VariableNames', varnames);
    speechelecs = [speechelecs; t2];
    clear t2
end
clear t2 s

df1 = 5 - 1; % degrees of freedom for numerator
df2 = length(Dvow.vowel) - 5; % degrees of freedom for denominator
fthresh = finv(1-0.0001, df1, df2);

tmp = linspace(fthresh, 30, numbins-2);
edges = [0 tmp 400];
for s=1:length(SIDs)
    SID = SIDs{s};
    if isfield(betaInfo, SID) && ~isempty(betaInfo.(SID).els)
        conds = discretize(fvals.(SID), edges); 
        conds(1, conds>1) = conds(1, conds>1)+1;
        
        % the second condition corresponds to speech responsive electrodes
        conds(1, intersect(find(conds==1),allidx_speech.(SID))) = 2;
        desel.(SID).elid = 1:length(fvals.(SID));
        desel.(SID).condition=conds;
    end        
end

% for single subject
desel.sz = [5, 25, repmat(55, 1, numbins)];
desel.cols = [ 0 0 0; 0 0 0 ; brewermap(numbins, 'YlGn')];

for s = SIDs    

    SID = s{1};
    [~] = plotNativeElec({SID}, desel, [datapath '/pt_data']);
    
    set(gcf,'Color','w');
    set(gcf, 'InvertHardcopy', 'off');
    title(SID);

    axes('Position',[.1 .1 .3 .3])
    p = pie([sum(desel.(SID).condition>2), sum(desel.(SID).condition==2)], [1 1]); 
    p(1).FaceColor = desel.cols(20, :);
    p(1).EdgeColor = 'none';
    p(3).FaceColor = [0.6 0.6 0.6];
    p(3).EdgeColor = 'none';
    p(2).FontSize=13; p(2).FontWeight = 'bold';
    p(4).FontSize=13; p(4).FontWeight = 'bold';

    %print(fullfile(['supp/' SID '_vowel.jpg']), '-djpeg', '-painters', '-r600')
end


clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections formant *model desel


%% Figure 1 - SUPP: Methods Pipeline

width = 450*2;
height = 2*200;
figure('Position', [100 800 width height]);
rc = [4 2];

% sound wave
subplot(rc(1), rc(2), 1);
corpus = 'dimex';
sentdet =getfield(load(['stim_info/out_sentence_details_' ...
    corpus '_all_loudness.mat']),'sentdet'); 
x = sentdet(4);

plot(x.sound, 'Color', [0.5 0.5 0.5]); set(gca, 'Visible', 'off');

% high gamma example
subplot(rc(1), rc(2), [2*rc(2)+1 3*rc(2)+1]);

filename = ['S2_B1_HilbAA_70to150_8band_all_0_mel_DIMEX_' ...
    'zflag_global_out_resp_log.mat'];
load([datapath '/pt_data/S2/dimex/' filename])

ctr=1;
for el=160:170    
    plot(ctr+out(1).resp(el, 1:350), 'k'); hold on;
    ctr=ctr+6;    
end
box off; yticks([]); %xticks([]);
xlabel('Time'); ylabel('Electrode');
set(gca, 'FontSize', 15); ylim([-7 70])
xlim([0 350]);

% plot spectrogram and other TRF features
% sentence onset + phonetic features
subplot(rc(1), rc(2), 2);
vertline(find(x.onsOff(1, :))', gcf, 'k', '-');
xlim([50 200]); yticks([]); xticks([]); box off;
ylabel({'Sentence', 'Onset'});
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',45,'VerticalAlignment','middle');
set(gca, 'FontSize', 13);

% peak rate
subplot(rc(1), rc(2), 2+rc(2))
plot(x.loudnessall(6, :), 'LineWidth', 3, 'Color', [0.5 0.5 0.5]);
xlim([50 200]); box off;
yticks([]); xticks([])
ylabel('peak rate');
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',45,'VerticalAlignment','middle');
set(gca, 'FontSize', 13);
 
% vowel onset
subplot(rc(1), rc(2), 2 + 2*rc(2));
vertline(x.vowelTimes(1, :)', gcf, 'k', '-');
xlim([50 200]); yticks([]); xticks([]); box off;
ylabel({'Vowel', 'Onset'});
set(gca, 'FontSize', 13);
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',45,'VerticalAlignment','middle');

% spectrogram
ax = subplot(rc(1), rc(2), 2 + 3*rc(2));
load('stim_info/mel_centerF.mat', 'binfrqs');
aud = binfrqs(1:end-2)+(binfrqs(3:end)-binfrqs(1:end-2))/2;

time = (1:size(x.aud, 2))/100;
imagesc(time, aud./1000, x.aud);
hold on;
set(gca, 'YDir', 'normal');
grays = colormap(ax, gray);
colormap(ax, flipud(gray));

idx =  x.formants(1, :)<2500 & x.formants(1, :)>0 & x.formants(2, :)<5000;
scatter(time(idx), x.formants(1, idx)./1000, 35, 'r', 'filled'); hold on;
scatter(time(idx), x.formants(2, idx)./1000, 35, 'r', 'filled'); hold on;

f0 = addF0(x, corpus, datapath);
idx = f0>0;
scatter(time(idx), f0(idx)./1000, 35,'b', 'filled');

yticks([1 4 8]);
ylim([0 8])
xlim([50/100 200/100]);
xticks(0.5:0.5:2);
xticklabels(split(num2str(0:0.5:1.5)))

ylabel({'Frequency', '(kHz)'})
hYLabel = get(gca,'YLabel');
set(hYLabel,'rotation',45,'VerticalAlignment','middle');

xlabel('Time (s)')
box off
set(gca, 'FontSize', 13);

% beta weights
figure; 
modelnames={'onset_maxDtL_maxDtLOnset_vowelOnset_F0_audF1F2'};
corpusStrf=loadMultModelStrf('S2', modelnames, corpus, ...
    [datapath '/pt_data'], 1, '');

strtfeat = 6;
f1bins = 15;
imagesc(corpusStrf{1}.meanStrf(strtfeat+f1bins:end, :, 168), ...
    'AlphaData', 0.8);
colormap(ax, flipud(pink))
ylabel({'Formant', 'Frequency Bins'});
title('Beta Weights');
xticks([1 60]); yticks([])
xticklabels({'0', '-0.6'});
xlabel('Time (s)');

cbh=colorbar;
cs=caxis;
set(cbh,'YTick',[cs(1) cs(2)], 'TickLabels', {'min', 'max'});
set(gca, 'FontSize', 13)

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* inflections ...
    formant *model desel

%% Figure 1 - SUPP: Unique Variance for F1 and F2

width = 300;
height = 300;
figure('Position', [100 800 width height]);

modelnames={'onset_maxDtL_maxDtLOnset_vowelOnset_F0_audF1', ...
    'onset_maxDtL_maxDtLOnset_vowelOnset_F0_audF2', ...
    'onset_maxDtL_maxDtLOnset_vowelOnset_F0_audF1F2', ...
    'onset_maxDtL_maxDtLOnset_vowelOnset_F0_aud', ...
    'onset_maxDtL_maxDtLOnset_vowelOnset_F0'};

load([datapath '/pt_data/out_elecs_voweltypeftest_bychan_anon.mat'], ...
    'allidx', 'fvals');
for s = SIDs
    SID = s{1};    
    corpusStrf=loadMultModelStrf(SID, modelnames, 'dimex', ...
        [datapath '/pt_data'], 1, '');

    idx = allidx.(SID);
    if ~any(cellfun(@(x) isempty(x), corpusStrf))
        univ_f1 = corpusStrf{3}.meanTestR.^2 -  corpusStrf{1}.meanTestR.^2;
        univ_f2 = corpusStrf{3}.meanTestR.^2 -  corpusStrf{2}.meanTestR.^2;
        %color = corpusStrf{3}.meanTestR.^2;
        color = fvals.(SID);
        disp(['added elecs: ' num2str(length(univ_f1(idx)))]);

        scatter(univ_f1(idx), univ_f2(idx), 55, color(idx), 'filled', ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'k'); 
        set(gca,'ColorScale','log');
        colormap(flipud(brewermap(64, 'Spectral')))
        hold on;
    else
        warning(['Did not include SID: ' SID]);
    end
end

xline(0, '-k');
yline(0, '-k');

ylabel('F2 bins unique variance');
ylim([-0.1 0.1])
yticks(-0.1:0.1:0.2)
xlabel('F1 bins unique variance');
xlim([-0.1 0.1])
xticks(-0.1:0.1:0.2)
cbh = colorbar;

ylabel(cbh, 'Vowel F-Stat');
caxis([5 60])
set(gca, 'FontSize', 13);

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* inflections ...
    formant *model desel


