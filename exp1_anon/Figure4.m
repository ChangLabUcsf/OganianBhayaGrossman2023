
% "STARTUP.M" NEEDS TO BE RUN BEFORE THE FOLLOWING CODE

%% --------------------- Speaker Normalization ----------------------------
infoStruct = {betaInfo, betaInfo_lte170, betaInfo_gt170}; 
%betaInfo_lt135, betaInfo_gt238};

inflections_f1 = loadF0Infl(1, infoStruct, SIDs);
inflections_f2 = loadF0Infl(2, infoStruct, SIDs);

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections* formant *model

%% Figure 4 (a, b) - Cross Speaker Formant Variation

figure('Position', [100 800 2*300 250]);

subplot(1, 2, 1);
% find medoids for separate speaker type
medoids_gt = plotVariance(Dvow.formantVals(:, Dvow.meanf0>170), ...
    Dvow.vowel(Dvow.meanf0>170), [], 0, 0);
medoids_lte = plotVariance(Dvow.formantVals(:, Dvow.meanf0<170), ...
    Dvow.vowel(Dvow.meanf0<170), [], 0, 0);

scatter(medoids_gt(:, 2), medoids_gt(:, 1), 105, ...
   getColors(1), 'filled'); hold on;
scatter(medoids_lte(:, 2), medoids_lte(:, 1), 105, ...
    getColors(1), 'filled', 's');

% figure style
set(gca, 'FontSize', 15);

set(gca, 'XDir', 'reverse',  'YDir', 'reverse');
xlabel('F2 (kHz)');  ylabel('F1 (kHz)');
xlim([500 3000]);  ylim([200 900]); 
xticks(1000:1000:3000); 
xticklabels(split(num2str(1:1:3)));

yticks(200:400:800);
yticklabels(split(num2str(0.2:0.4:0.8)));
legend('off');

% formant variation
subplot(1, 2, 2);
medoids = plotVariance(Dvow.formantVals,  Dvow.vowel, [], 0, 0);
med = cell2mat(arrayfun(@(x) medoids(x, 1:2), ...
    Dvow.vowelType, 'UniformOutput', false)');
med_lte = cell2mat(arrayfun(@(x) medoids_lte(x, 1:2), ...
    Dvow.vowelType, 'UniformOutput', false)');
med_gt = cell2mat(arrayfun(@(x) medoids_gt(x, 1:2), ...
    Dvow.vowelType, 'UniformOutput', false)');

var = arrayfun(@(x) pdist([Dvow.formantVals(1:2, x)'; med(x, :)]), ...
    1:size(Dvow.formantVals, 2));

lte_var = arrayfun(@(x) pdist([Dvow.formantVals(1:2, x)'; med_lte(x, :)]), ...
    find(Dvow.meanf0<170)); 
gt_var = arrayfun(@(x) pdist([Dvow.formantVals(1:2, x)'; med_gt(x, :)]), ...
    find(Dvow.meanf0>170)); 

bw = 0.05;
violin(rmoutliers(var./1000)', 'x', 1:2, 'facecolor', [0.1 0.1 0.1], ...
    'edgecolor','none', 'bw', bw, 'medc', []); hold on;
violin(rmoutliers(lte_var./1000)', 'x', 2:3, 'facecolor', [0.4 0.4 0.4], ...
    'edgecolor','none', 'bw', bw, 'medc', []); hold on;
violin(rmoutliers(gt_var./1000)', 'x', 3:4, 'facecolor', [0.8 0.8 0.8], ...
    'edgecolor','none', 'bw', bw, 'medc', []); hold on;
xticks(1:3);
xticklabels({'all', 'F0<170', 'F0>170'})

% legend({'All', 'Separated by F0'})
xlim([0 3.5]);
ylim([-0.1 0.9]);
set(gca, 'FontSize', 15);
ylabel({'Distance', 'to Medoid (kHz)'});
box off;

[~, p] = ttest2(var, lte_var);
disp(['ttest between variance vs. F0<170Hz variance: ' num2str(p, 3)]);
if p<0.001
    sigline([1, 2],'p<.001',0.6);
end

[~, p] = ttest2(var, gt_var);
disp(['ttest between variance vs. F0>170Hz variance: ' num2str(p, 3)]);
if p<0.001
    sigline([1,3],'p<.001',0.7);
end
legend off


clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections* formant *model

%% Figure 4 (d, e) - Speaker Normalization Example Electrode 1 

elecs = containers.Map;

% figure electrodes
elecs('S6') = 33; % 83
elecs('S1') = 71; % 87

corpus='dimex';

numel = length(cell2mat(values(elecs)));
figure('Position', [100 800 300*4 2*250]);

subplt = 1;
for key = keys(elecs)    
    SID = key{1};
    Dvow = addtoDD(Dvow, corpus, bef, aft, {SID});
       
    for el = elecs(SID)
        
        tmpidx = find(betaInfo.(SID).els==el);

        bins = [13 13];
        ax = subplot(numel, 3, subplt);
        plotImagescResponse(Dvow, SID, el,  1, ax, getColors(1));    
        
        % beta fitting
        cols = {'b', 'r'};
        ax=nan(4, 1);
                
        for f=1:2            
            % sigmoid definition
            fsigm = @(param,xval) ...
                       param(1)+(param(2)-param(1))./...
                       (1+10.^((param(3)-xval)*param(4)));
            
            linestyle = {'-', '--'};
            alpha = [0.8, 0.4];
            ctr = 1; 
            for i = {betaInfo_gt170, betaInfo_lte170} %{betaInfo_gt170, betaInfo_lte170} % 
                % plot sigmoid fit and scatter the real betas
                info = i{1};
                
                ax(f) = subplot(numel, 3, subplt + f);
                x = info.x(f, :); 
               
                y_norm = normalize(info.(SID).y(tmpidx, f, :));
                scatter(betaInfo.x(f, :), squeeze(y_norm), 25, 'k', ...
                    'filled', 'DisplayName', ['F' num2str(f)], ...
                    'MarkerFaceAlpha', alpha(ctr)); hold on;

                slope = info.(SID).beta_lin(tmpidx, f, 1);
                
                y_sigm = normalize(fsigm(info.(SID).beta_sigm(tmpidx, f, :), x));                
                h=plot(x, squeeze(y_sigm), [cols{(slope<0)+1} linestyle{ctr}], ...
                    'HandleVisibility', 'off', 'LineWidth', 2);

                infl = info.(SID).beta_sigm(tmpidx, f, 3);                
                infl_y = interp1(x(~isnan(x)), y_sigm(~isnan(x)), infl);
                scatter(infl, infl_y, 135, 'k', 'x');
                xlim([min(x) max(x)]);
                ctr = ctr+1;       
                
                clear y_sigm ys  
            end
            ylim([-2 2]);
            yticks(-2:2:2);
            
            % remove scientific notation
            ax = ancestor(h, 'axes');
                                        
            xs = xticks;
            xticklabels(split(num2str(xs./1000)));
            xlabel(['F' num2str(f) ' (kHz)']);
            if f == 1               
                ylabel({'Norm. Beta', 'Weights (a.u.)'});
            end
            
            title({['F0>170 R^2: ' num2str(betaInfo_gt170.(SID).rsq(1, tmpidx, f), 2)], ...
                ['F0<170 R^2: ' num2str(betaInfo_lte170.(SID).rsq(3, tmpidx, f), 2)]}, ...
                'FontWeight', 'normal');
            set(gca, 'FontSize', 15);           
        end
        subplt = subplt + 3;    
    end
end

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    *model inflections* formant numel subplt

%% Figure 4 (f, h, g) - Speaker Normalization Inflection Point Shifts

fig = figure('Position', [100 800 300*3 1*250], 'Renderer', 'Painters');
colormap_rsq = brewermap(5, 'YlGn');

% average acoustic difference between two speaker groups
acs_diff = mean(Dvow.formantVals(:, Dvow.meanf0>170), 2) - ...
    mean(Dvow.formantVals(:, Dvow.meanf0<170), 2);

rsqth = 0.1;
infls = {inflections_f1, inflections_f2};
lims = [0.2 0.75; 0.7 2.5];
for formant = 1:2
    ax = subplot(1, 2, formant);
    switch formant
        case 1
            range = 0.2:0.2:0.8;
        case 2
            range = 1:0.5:2.5;
    end
    
    inflections = infls{formant};    
    idx = inflections.sr_lte170>rsqth & inflections.sr_gt170>rsqth & ...
        inrange(inflections.infl_lte170./1000, lims(formant, :)) & ...
        inrange(inflections.infl_gt170./1000, lims(formant, :));
    
    % 2D colormap
    numColors = 21;
    [cm2d] = ecog_2dcm(round(numColors/2), [1 1 1; 0 1 1 ; 1 1 1; 1 1 0]);
      
    edges = linspace(0, 1, numColors);
    [col1, ~] = discretize(inflections.sr_lte170(idx), edges);
    [col2, ~] = discretize(inflections.sr_gt170(idx), edges);
    color = cell2mat(arrayfun(@(i) squeeze(cm2d(col1(i), col2(i), :)), ...
        1:sum(idx, 'omitnan'), 'UniformOutput', false))';    
        
    sz = 65; 
    xticks(-10:10:10);
    xticklabels(split(num2str(edges(1:10:25))));
    yticks(-10:10:10);
    yticklabels(split(num2str(edges(1:10:25))));
    set(gca, 'FontSize', 13);
    ylabel('F0 > 170 Hz R^2');
    xlabel('F0 < 170 Hz R^2');
    
    set(0, 'CurrentFigure', fig);
    z = inflections.el(idx);%SID(idx);
    scatter3(inflections.infl_gt170(idx)./1000, inflections.infl_lte170(idx)./1000, z, ...
        sz, color, 'filled', 'MarkerFaceAlpha', 0.9, 'MarkerEdgeColor', 'k'); 
    
    colormap(ax, colormap_rsq)
    labels = {'>170', '<=170'};
    xlabel([labels{1} ' Inflection (kHz)']); 
    ylabel([labels{2} ' Inflection (kHz)']);

    % figure style
    view(2); legend('show');
    legend('off');
    
    yticks(range);
    xticks(range);
    ylim(lims(formant, :));
    xlim(lims(formant, :));

    % reference lines
    h=refline(1, 0);
    h.Color = 'k';
    h.LineWidth = 2;   
    h.HandleVisibility='off'; 
    h=refline(1, -acs_diff(formant)/1000);
    h.Color = 'k';
    h.LineStyle = '--'; 
    h.LineWidth = 1.5;   
    h.HandleVisibility='off'; 
    set(gca, 'FontSize', 15);
    ylim(lims(formant, :));
    xlim(lims(formant, :));
    grid off
    title(['n = ' num2str(sum(idx, 'omitnan'))]);
end

total_elec = (infls{1}.sr_lte170>rsqth & infls{1}.sr_gt170>rsqth & ...
        inrange(infls{1}.infl_lte170./1000, lims(1, :)) & ...
        inrange(infls{1}.infl_gt170./1000, lims(1, :))) | ...
        (infls{2}.sr_lte170>rsqth & infls{2}.sr_gt170>rsqth & ...
        inrange(infls{2}.infl_lte170./1000, lims(2, :)) & ...
        inrange(infls{2}.infl_gt170./1000, lims(2, :)));
disp(['Total electrodes in F & G: ' num2str(sum(total_elec)) '' ...
    ', percentage: ' num2str((sum(total_elec)./height(infls{1}))*100) '%']);


% Figure 2 - Speaker Normalization over F1 and F2 (Co-occurrence)
infl_f1 = inflections_f1;
infl_f2 = inflections_f2;

% make sure both rsq meet threshold
sr_f1 = infl_f1.('sr_lte170') > rsqth & infl_f1.('sr_gt170') > rsqth;
sr_f2 = infl_f2.('sr_lte170') > rsqth & infl_f2.('sr_gt170') > rsqth;
sr = sr_f1 & sr_f2;

% find inflection difference
infldiff_f1 = infl_f1.('infl_gt170') - infl_f1.('infl_lte170');
infldiff_f2 = infl_f2.('infl_gt170') - infl_f2.('infl_lte170');

fig2 = figure('Position', [100 800 300*1 1*250], 'Renderer', 'Painters');
set(0, 'CurrentFigure', fig2);
ax = subplot(1, 1, 1);
z = infl_f1.el; 

numColors = 21;
[cm2d] = ecog_2dcm(round(numColors/2), [1 1 1; 1 0 0 ; 1 1 1; 0 0 1]);
edges = linspace(0, 1, numColors);

[col1, edges_lte] = discretize(max([infl_f1.('sr_gt170')(sr) ...
    infl_f1.('sr_lte170')(sr)], [], 2, 'omitnan'), edges);
[col2, edges_gt] = discretize(max([infl_f2.('sr_gt170')(sr) ...
    infl_f2.('sr_lte170')(sr)], [], 2, 'omitnan'), edges);
color = cell2mat(arrayfun(@(i) squeeze(cm2d(col1(i), col2(i), :)), ...
    1:sum(sr, 'omitnan'), 'UniformOutput', false))'; 

sz = 75; 
xticks(-10:10:10);
xticklabels(split(num2str(edges(1:10:25))));
yticks(-10:10:10);
yticklabels(split(num2str(edges(1:10:25))));
set(gca, 'FontSize', 13);
ylabel('F2 max R^2');
xlabel('F1 max R^2');

set(0, 'CurrentFigure', fig2);

% within range
idx = inrange(infldiff_f1./1000, [-2 2]) & inrange(infldiff_f2./1000, [-4 4]);
scatter3(infldiff_f1(sr)./1000, infldiff_f2(sr)./1000, z(sr), 65, color, ... %  sz(sr)
    'filled', 'MarkerFaceAlpha', 0.8); % 'MarkerEdgeColor', 'k'

xlim([-0.2 0.4]);
ylim([-1.5 2.5]);
yline(0, '-k', 'LineWidth', 1.2);
xline(0, '-k', 'LineWidth', 1.2);
yline(acs_diff(2)/1000, '--k', 'LineWidth', 1.5);
xline(acs_diff(1)/1000, '--k', 'LineWidth', 1.5);
xticks(-0.2:0.2:0.4);
yticks(0:2:2);
view(2);
grid off;

set(gca, 'FontSize', 15);
ylabel({'Infl. Difference', '(F2, kHz)'});
xlabel({'Infl. Difference', '(F1, kHz)'});

title(['n = ' num2str(sum(sr&idx, 'omitnan'))]);

% Running the statistical models
ctr = 1;
lim = [200 1000; 200 3000];
for j = {infl_f1, infl_f2}
    i = j{1};
    numel = height(i);
    % el = cellstr(num2str(repmat(inflections_f1.el, 3, 1)));
    el = repmat(i.el, 2, 1);
    infl = [i.infl_lte170; i.infl_gt170];
    pitch = [-1*ones(numel, 1); ones(numel, 1);];
    rsq = [i.sr_lte170; i.sr_gt170];
    sid = repmat(i.SID, 2, 1);

    % add speaker normalization regression
    tbl = table(infl, pitch, rsq, el, sid, ...
        'VariableNames',{'Infl','Pitch','Rsq','El', 'Subj'});
    
    % within 
    idx = tbl.Infl<lim(ctr, 1) | tbl.Infl>lim(ctr, 2);
    % use original rsqthresh (0.1 for consistency)
    % significance of mixed-effect is not altered
    idx = idx | tbl.Rsq<0.15; 
    tbl(idx, :) = [];
    
    tbl.Rsq = zscore(tbl.Rsq); 
    
    lme = fitlme(tbl,'Infl~Rsq*Pitch+(1|Subj)+ (1|El:Subj)');
    disp(['-------------------------F' num2str(ctr) '------------------']);
    disp(lme);
    ctr = ctr+1;
end

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections* formant *model

%% Figure 4 - SUPP: Marginals for F1 and F2, across and within speaker types

width = 455*2;
height = 200*2;
figure('Position', [100 800 width height]);

bm = getColors(1);

lims = {[100 900], [100 3000]};

pitchlabels = {'all', 'lt170', 'gt170'};
for pitch = 1:3 % 
    switch pitch
        case 1
           trials = Dvow.meanf0 > 0;
        case 2
           trials = Dvow.meanf0 < 170;
        case 3
           trials = Dvow.meanf0 > 170;
    end
    for formant=1:2
        subplot(2, 3, (formant-1)*3+pitch);
        for v = 1:5
            [d1_ac, EDGES] = discretize(Dvow.formantVals(formant, ...
                Dvow.vowelType == v & trials), 30);
            plot(EDGES(1:end-1)+diff(EDGES)./2, histcounts(d1_ac), ...
                'LineWidth', 1.5, 'Color', bm(v, :)); hold on;        
        end
        xlim(lims{formant});
        xlabel(['Formant ' num2str(formant) 'Hz']);
        if formant<2, title(['pitch group: ' pitchlabels{pitch}]); end
    end
end
legend(unique(Dvow.vowel));

% alternatively, show all five vowels with different pitch marginals
width = 455*2;
height = 200*2;
figure('Position', [100 800 width height]);

lims = {[200 900], [600 3000]};
ls = {'-', '--', ':'};
for v = 1:5    
    for formant=1:2
        subplot(2, 5, (formant-1)*5+v);
        for pitch = 1:3 % 'all', 'lt170', 'gt170'
            switch pitch
                case 1
                   trials = Dvow.meanf0 > 0;
                case 2
                   trials = Dvow.meanf0 < 170;
                case 3
                   trials = Dvow.meanf0 > 170;
            end

            [d1_ac, EDGES] = discretize(Dvow.formantVals(formant, ...
                Dvow.vowelType == v & trials), 30);
            plot(EDGES(1:end-1)+diff(EDGES)./2, histcounts(d1_ac), ...
                'LineWidth', 1.5, 'Color', bm(v, :), 'LineStyle', ls{pitch}); hold on;        
        end
        xlim(lims{formant});
        xlabel(['Formant ' num2str(formant) 'Hz']);
    end
end
legend({'all', '<170', '>170'});

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* inflections*...
    formant *model


%% Figure 4 - SUPP: Latency of Normalizing vs. Non-normalizing Electrodes

infls = {inflections_f1, inflections_f2};

column = {'f1', 'f2'};
hFig=figure('Renderer', 'painters', 'Position', [100 800 350*2 250]); 
for formant = 1:2 
    inflections = infls{formant};
    for s = SIDs
        SID = s{1};

        corpStrf = loadMultModelStrf(SID, {model}, 'dimex', ...
            [datapath '/pt_data'], 1);
        corpStrf = corpStrf{1};

        % find the beginning of the aud feature, assuming all other features
        % are 1-dimensional
        strtfeat = find(cellfun(@(x) strcmp(x, 'aud'), corpStrf.featureNames));
        [~, maxidx] = max(squeeze(mean(corpStrf.meanStrf(strtfeat:end, :, :), ...
            'omitnan')), [], 'omitnan');

        % only rows where SID is same as current SID
        subtable = inflections(inflections.SID == str2double(SID(2:end)), :);

        infl_lte = subtable.('infl_lte170');
        infl_gt = subtable.('infl_gt170');  
        sr = subtable.('sr_gt170') > 0 & subtable.('sr_lte170') > 0;
        elecs = subtable(sr, 'el');
        latency = maxidx(elecs.el);
        infl_diff = infl_gt(sr) - infl_lte(sr);
        
        set(0,'CurrentFigure',hFig)
        subplot(1, 2, formant);
        color = max([subtable.('sr_gt170') subtable.('sr_lte170')], [], 2, 'omitnan');
        scatter(infl_diff, latency*0.01, 45, color(sr), 'filled', ...
            'MarkerEdgeColor', 'k'); hold on;        
    end   
end

% formatting
for formant = 1:2
    subplot(1, 2, formant);
    
    if formant == 1
        xlim([-500 500]);
        xticks(-500:500:500)
    else
        xlim([-2000 2000]);
        xticks(-2000:1000:2000)
    end
    ylim([0 0.25]);
    yticks(0:0.1:0.3)

    colormap(brewermap(30, 'PRGn'));
    cbh=colorbar;
    caxis([0 1]);
    set(cbh,'YTick',0:0.5:1);
    ylabel(cbh, 'Sigm R^2');

    ylabel('Latency of Peak Response (s)');
    if i == 1
        xlabel('Inflection Difference');
    end
    set(gca, 'FontSize', 15);

    h=xline(0);
    h.Color = 'k';
end

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections* formant *model

%% Figure 4 - SUPP: Inflection Location as Speaker Mean Pitch Increases 

% order is infl, infl_lte, infl_gt
infoStruct = {betaInfo_lt135, betaInfo_bt135_200, betaInfo_bt200_236}; %, betaInfo_lte170, betaInfo_gt170,
pitchStr = {'<135 Hz', '135-200 Hz', '200-236 Hz'};
xlin = [135, 200, 236];
%infoStruct = {betaInfo_lt135, betaInfo_bt135_238, betaInfo_gt170};

thresh = 0.05;
inflections_f1 = loadF0Infl(1, infoStruct, SIDs);
inflections_f2 = loadF0Infl(2, infoStruct, SIDs);
infls = {inflections_f1, inflections_f2};

figure;
h1 = subplot(2, 2, 1);
h2 = subplot(2, 2, 2);
sbs = [h1, h2];

% scatter of inflection points per group
for f = 1:2 
    set(gcf,'CurrentAxes', sbs(f));
    
    % boxplot
    infl = infls{f};
    idx = infl.sr>thresh & infl.sr_lte170>thresh & ...
        infl.sr_gt170>thresh;

    id1 = infl.sr>thresh; % significant in first speaker subset
    id2 = infl.sr_lte170>thresh; % significant in second speaker subset
    id3 = infl.sr_gt170>thresh; % significant in third speaker subset
    disp(['---------------------- F' num2str(f) ' ----------------------']);
    disp(['sigmoid rsq (lt, btb, gt) over 0.5: ' num2str(sum(id1, 'omitnan')) ', ' ...
        num2str(sum(id2, 'omitnan')) ', ' num2str(sum(id3, 'omitnan'))]);
    disp(['intersect: ' num2str(sum(all([id1, id2, id3], 2), 'omitnan'))]);
    disp(['numel: ' num2str(sum(idx, 'omitnan'))])
    
    % Not useful, better to use pitch as different stepwise predictors (1, 2, 3)
%     disp(['-------------------- F' num2str(f) ' ftest --------------------']);
%     [~,p1] = ttest2(infl.infl(idx), infl.infl_lte170(idx));
%     [~,p2] = ttest2(infl.infl_lte170(idx), infl.infl_gt170(idx));
%     [~,p3] = ttest2(infl.infl(idx), infl.infl_gt170(idx));
%     disp(['pair (1, 2): ' num2str(p1)])
%     disp(['pair (2, 3): ' num2str(p2)])
%     disp(['pair (1, 3): ' num2str(p3)])    

    boxplot(sbs(f), [infl.infl(idx) infl.infl_lte170(idx) infl.infl_gt170(idx)], ...
        'symbol','');
    set(findobj(gca,'type','line'),'linew',0.8, 'color', 'r'); hold on;
    
    %  overlaid on scatter
    for i = 1:height(inflections_f1)  
        infl = infls{f};
        if idx(i) % only plot if part of barplot
            y = [infl.infl(i) infl.infl_lte170(i) infl.infl_gt170(i)];
            alpha = [infl.sr(i) infl.sr_lte170(i) infl.sr_gt170(i)];
            alpha(alpha<0 | isnan(alpha)) = 0;
            r = nan(3, 1);
            for x = 1:3   
                r(x) = rand;
                scatter(x-0.1+r(x)*0.2, y(x), 35, 'k', 'filled', 'MarkerFaceAlpha', ...
                    alpha(x)); hold on;        
            end
    %         if min(alpha)>0.8
    %             lh = plot(h,[1:3]'-0.1+r*0.2, y, 'LineWidth', min(alpha)); hold on;
    %             lh.Color=[0,0,0,2*min(alpha)-1];
    %         end    
        end
    end
    
    % figure style    
    xlim([0.6 3.4]);
    xticks(1:3);   
    xticklabels(pitchStr);
    ylabel(['Infl F' num2str(f) ' (Hz)']);
    set(gca, 'FontSize', 15);
    switch f
        case 1
            yticks(200:200:1000);
            ylim([300 700]);
        case 2
            yticks(1000:1000:3000)
            ylim([800 3000]);            
    end
    box off;
end

subplot(2, 2, [3 4]);
histogram(Dvow.meanf0, 'EdgeColor', 'none', 'FaceAlpha', 0.5, ...
    'FaceColor', [0.5 0.5 0.5]);
xline(xlin, 'LineWidth', 3, 'Color', 'r');

% pairwise difference
% figure;
% h1 = subplot(1, 2, 1);
% h2 = subplot(1, 2, 2);
% sbs = [h1, h2];
% for f = 1:2
%     set(gcf,'CurrentAxes', sbs(f));
%     % boxplot overlaid on scatter
%     infl = infls{f};
%     idx = infl.sr>thresh & infl.sr_lte170>thresh & ...
%         infl.sr_gt170>thresh;
%     
%     y = [inflections_f1.infl_lte170(idx) - inflections_f1.infl(idx) ...
%     inflections_f1.infl_gt170(idx) - inflections_f1.infl_lte170(idx)];
%     boxplot(y);
%     
%     set(findobj(gca,'type','line'),'linew',0.8, 'color', 'r'); hold on;    
%     yline(0, 'LineStyle', '--', 'Color', 'k');
%     
%     switch f
%         case 1
%             ylim([-300 300]);
%         case 2    
%             ylim([-1500 1500]);
%     end    
% end

% single electrode example

% figure electrodes
elecs = containers.Map;
%elecs('EC214') = 3;

% example electrpde
elecs('EC203') = []; %129
%elecs('EC214') = 69;

corpus='dimex';
addpath othercolor

numel = length(cell2mat(values(elecs)));
formant = 2;

subplt = 1;
for key = keys(elecs)    
    SID = key{1};
       
    for el = elecs(SID)   
        figure('Position', [100 800 300*4 2*250]);
        tmpidx = find(betaInfo.(SID).els==el);
        
        % beta fitting
        cols = {'b', 'r'};
        ax=nan(4, 1);
            
        % sigmoid definition
        fsigm = @(param,xval) ...
                   param(1)+(param(2)-param(1))./...
                   (1+10.^((param(3)-xval)*param(4)));

        labels = pitchStr;
        linestyle = {'-', '--', '-.'};
        alpha = [0.9, 0.6, 0.3];
        ctr = 1; 
        for i = infoStruct
            % plot sigmoid fit and scatter the real betas
            info = i{1};
            x = info.x(formant, :); 

            y_norm = normalize(info.(SID).y(tmpidx, formant, :));
            scatter(betaInfo.x(formant, :), squeeze(y_norm), 25, 'k', ...
                'filled', 'MarkerFaceAlpha', alpha(ctr), ...
                 'HandleVisibility', 'off'); hold on;

            slope = info.(SID).beta_lin(tmpidx, formant, 1);

            y_sigm = normalize(fsigm(info.(SID).beta_sigm(tmpidx, formant, :), x));                
            h=plot(x, y_sigm, [cols{(slope<0)+1} linestyle{ctr}], ...
                'LineWidth', 2, 'DisplayName', labels{ctr});

            infl = info.(SID).beta_sigm(tmpidx, formant, 3);                
            infl_y = interp1(x(~isnan(x)), y_sigm(~isnan(x)), infl);
            scatter(infl, infl_y, 135, 'k', 'x', 'HandleVisibility', 'off');
            xlim([min(x) max(x)]);
            ctr = ctr+1;       

            clear y_sigm ys  
        end
        ylim([-2 2]);
        yticks(-2:2:2);

        xs = xticks;
        xticklabels(split(num2str(xs./1000)));
        xlabel(['F' num2str(formant) ' (kHz)']);              
        ylabel({'Norm. Beta', 'Weights (a.u.)'});
      
        set(gca, 'FontSize', 15);           
        subplt = subplt + 3; 
    end
end

ctr = 1;
rng = [100 1000; 100 4000];
for j = infls
    i = j{1};
    numel = height(i);
    % el = cellstr(num2str(repmat(inflections_f1.el, 3, 1)));
    el = repmat(i.el, 3, 1);
    infl = [i.infl; i.infl_lte170; i.infl_gt170];
    pitch = [-1*ones(numel, 1); 0*ones(numel, 1); ones(numel, 1)];
    rsq = [i.sr; i.sr_lte170; i.sr_gt170];
    % add speaker normalization regression
    tbl = table(infl, pitch, rsq, el, 'VariableNames',{'Infl','Pitch','Rsq','Elec'});
    
    idx = tbl.Infl<rng(ctr, 1) | tbl.Infl>rng(ctr, 2);
    idx = idx | tbl.Rsq<0.15 |isnan(tbl.Infl); 
    tbl(idx, :) = [];
    
    %tbl.Pitch = zscore(tbl.Pitch);
    tbl.Rsq = zscore(tbl.Rsq);
    
    lme = fitlme(tbl,'Infl~Rsq*Pitch+(1+Pitch|Elec)');
    disp(['-------------------------F' num2str(ctr) '---------------------------']);
    disp(lme);
    ctr = ctr+1;
end

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections* formant *model


%% -------------------------- Functions -----------------------------------
function [inflections] = loadF0Infl(formant, infoStruct, SIDs)
   
    % set properties of analysis
    % variable so we can use different speaker betas
    % infoStruct = {betaInfo, betaInfo_lte170, betaInfo_gt170}; 
    
    inflections =  array2table(zeros(0,8), 'VariableNames', ...
        {'SID', 'el', 'sr', 'infl', 'sr_lte170', 'infl_lte170', 'sr_gt170', 'infl_gt170'});

    for s = 1:length(SIDs)
        SID=SIDs{s};
        ctr = 1;

        if isfield(infoStruct{1}, SID) && ~isempty(infoStruct{1}.(SID).els)   
            sr = nan(length(infoStruct), length(infoStruct{1}.(SID).els))';
            infl = nan(length(infoStruct), length(infoStruct{1}.(SID).els))';
            for i = infoStruct
                info = i{1};
                if isfield(info, SID) && ~isempty(info.(SID).els)           
                % difference between sigmoid and linear rsq (for both formants)   

                    % significant linear vs significant rsq
                    sr(:, ctr) = squeeze(info.(SID).rsq(3, :, formant))'; 

                    % inflection points
                    infl(:, ctr)=info.(SID).beta_sigm(:, formant, 3);                     
                    clear sigrsq* linrsq* 
                end
                ctr = ctr + 1;
            end
            sids = repmat(str2double(SID(2:end)), length(infoStruct{1}.(SID).els), 1);

            if ~isempty(infoStruct{1}.(SID).els)
                t2 = array2table([sids infoStruct{1}.(SID).els sr(:, 1) infl(:, 1) ...
                            sr(:, 2) infl(:, 2) sr(:, 3) infl(:, 3)], ...
                            'VariableNames',  {'SID', 'el', 'sr', 'infl', 'sr_lte170', 'infl_lte170', ...
                            'sr_gt170', 'infl_gt170'});
                inflections = [inflections; t2];
            end
        end
    end
end