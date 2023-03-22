
% "STARTUP.M" NEEDS TO BE RUN BEFORE THE FOLLOWING CODE

%% Figure 2 - Load inflection data

[inflections]=loadInfl(betaInfo, SIDs);

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* inflections ...
    formant *model

%% Figure 2 (b, c, d, e) - Example RFs

%Linear: EC105 - 168
%Sigmoid: EC214 - 54
elecs = containers.Map;
elecs('S1') = 71; %105 - 152, 168;
elecs('S7') = 54;
corpus='dimex';

numel = length(cell2mat(values(elecs)));
fig = figure('Position', [100 800 300*3 1 + numel*250]);

ctr = 1;
for key = keys(elecs)    
    SID = key{1};
    Dvow = addtoDD(Dvow, corpus, bef, aft, {SID});
       
    for el = elecs(SID)

        % high gamma vowel map 
        ax = subplot(numel, 3, ctr);
        plotImagescResponse(Dvow, SID, el,  1, ax, getColors(1));    
        
        tmpidx = find(betaInfo.(SID).els==el);        
        
        % beta fitting
        cols = {'b', 'r'};
        ax=nan(4, 1);
                
        for f=1:2       
            ax(f) = subplot(numel, 3, ctr+f);
            x = betaInfo.x(f, :); 
                
            y_norm = normalize(betaInfo.(SID).y(tmpidx, f, :));
            scatter(betaInfo.x(f, :), squeeze(y_norm), 25, 'k', ...
                'filled', 'DisplayName', ['F' num2str(f)], ...
                'MarkerFaceAlpha', 0.6); hold on;
            
            f1 = normalize(polyval(squeeze(betaInfo.(SID).beta_lin(tmpidx, f, :)), x));          

            slope = betaInfo.(SID).beta_lin(tmpidx, f, 1);
            plot(x,f1, [cols{(slope<0)+1} '-'], 'HandleVisibility', 'off', ...
                'LineWidth', 1); hold on;                  

            % plot sigmoid fit and scatter the real betas
            fsigm = @(param,xval) ...
                   param(1)+(param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)));

            y_sigm = normalize(fsigm(betaInfo.(SID).beta_sigm(tmpidx, f, :), x));                        
            plot(x,y_sigm, [cols{(slope<0)+1} '--'], 'HandleVisibility', ...
                'off', 'LineWidth', 2);

            infl = betaInfo.(SID).beta_sigm(tmpidx, f, 3);
            infl_y = interp1(x(~isnan(x)), y_sigm(~isnan(x)), infl);
            scatter(infl, infl_y, 35, cols{(slope<0)+1}, 'filled');
            xlim([min(x) max(x)]);

            ylim([-2 2]) 
            yticks(-2:2:2)
            clear y_sigm ys  
            
            % remove scientific notation
            linkaxes(ax, 'y');

            if f == 1
                ylabel('Beta Weights (a.u.)')
            end
            
            xs = xticks;
            xticklabels(split(num2str(xs./1000)));
            xlabel(['F' num2str(f) ' (kHz)']);
            
            linrsq = betaInfo.(SID).rsq(1, tmpidx, f);
            sigrsq = betaInfo.(SID).rsq(3, tmpidx, f);
            
            if sigrsq - linrsq > 0.05                
                title({['\rm Linear Fit R^2: ' num2str(linrsq, 2)], ...
                    ['\bf Sigmoid Fit R^2: ' num2str(sigrsq, 2)]});
            else
                title({['\bf Linear Fit R^2: ' num2str(linrsq, 2) ], ...
                    ['\rm Sigmoid Fit R^2: ' num2str(sigrsq, 2)]});
            end
            set(gca, 'FontSize', 13);
        end    
        ctr = ctr + 3;        
    end
end

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections formant *model

%% Figure 2 (f, g) - Sigmoidal vs. Linear Fit Over All Electrodes

figure('Renderer', 'painters', 'Position', [100 800 300*2 2*250]); 
    
column = {'f1', 'f2'};
rc = [2, 2];

for f = 1:2
    subplot(rc(1), rc(2), f*2);
    sr = inflections.(['sr_' column{f}]);
    lr = inflections.(['lr_' column{f}]);
    slope =  inflections.(['sl_' column{f}]);
    
    z = inflections.el;

    pos = sr > 0 | lr > 0;
    sig = sum(sr > 0 & sr > lr, 'omitnan');
    lin = sum(lr > 0 & lr > sr, 'omitnan');
    rdiff = abs(sr(pos)-lr(pos));

    % red for negative, blue for positive
    RB = [1 0 0; 0 0 1]; 
    colors = cell2mat(arrayfun(@(x) RB(x+1, :), slope(pos), ...
        'UniformOutput', false));                

    % scatter of sigmoidal and linear rsquared values
    scatter3(lr(pos), sr(pos), z(pos), 25, colors, 'filled', ...
            'MarkerFaceAlpha', 0.5, ...
            'HandleVisibility', 'off'); hold on;    
        
    disp(['---------------' column{f} '---------------']); 
    disp(['Plotting ' num2str(sum(pos, 'omitnan')) '/' num2str(length(pos)) ...
        ' electrodes']);    
    disp(['Distr of sig and lin ' num2str(sig) ' - ' num2str(lin)]);
    disp(['Percent sig ' num2str(100*(sig/(lin+sig)))]);
end       

% setting figure style
for i=1:2
    subplot(rc(1), rc(2), i*2);
   
    xlabel('linear R^2'); ylabel('sigmoid R^2');
    
    % reference lines
    xlim([-0.5 1]); 
    ylim([-0.5 1]);
    h=refline(1, 0); 
    h.LineStyle='--'; 
    h.Color='k';
    h.HandleVisibility='off';
    yline(0, 'Color', 'k', 'HandleVisibility', 'off'); 
    xline(0, 'Color', 'k', 'HandleVisibility', 'off');
    
    set(gca, 'FontSize', 13); 
    view(2); 
    grid off;
end

% Figure 1 - Monotonic histogram
maxidx = nan(height(inflections), 2);
minidx = nan(height(inflections), 2);
ctr = 1;
for s = SIDs
    SID = s{1};
    % el x formant x elem   
    if isfield(betaInfo.(SID), 'modelrsq')
        rsq = betaInfo.(SID).modelrsq;
        [~, idx] = max(betaInfo.(SID).y(rsq>0.02, 1:2, :), [], 3);
        maxidx(ctr:ctr+size(idx, 1)-1, :) = idx;

        % aggregating slope, electrode, and subject information
        sl(ctr:ctr+size(idx, 1)-1, :) = betaInfo.(SID).beta_lin(rsq>0.02, 1:2, 1)>0;
        el(ctr:ctr+size(idx, 1)-1, :) = betaInfo.(SID).els(rsq>0.02);
        subj(ctr:ctr+size(idx, 1)-1) = repmat(str2double(SID(2:end)), 1, size(idx, 1));

        [~, idx] = min(betaInfo.(SID).y(rsq>0.02, 1:2, :), [], 3);
        minidx(ctr:ctr+size(idx, 1)-1, :) = idx;

        ctr = ctr + length(idx);
    end
end

xlims = [0.2 0.8; 0.6 3];
fb = [{0.0:0.0488:1} {0.0:0.3:3.5}]; % bandwith
b = [0.06, 0.17];
for f = 1:2
    subplot(rc(1), rc(2), 2*f-1);
    
    % positive slope
    id = maxidx(sl(:, f), f);
    
    pos_subj = subj(sl(:, f));
    pos_subj(isnan(id)) = [];

    pos_el = el(sl(:, f));
    pos_el(isnan(id)) = [];

    id(isnan(id))=[];

    pos_slope = betaInfo.x(f, id)./1000;    
    histogram(pos_slope, fb{f}, 'FaceColor', 'b', 'FaceAlpha', 0.3, ...
        'EdgeColor', 'none', 'Normalization', 'probability'); hold on;
    
    yyaxis right
    [y_p, x_p] = ksdensity(pos_slope, 'BandWidth', b(f));
    plot(x_p, y_p, 'LineWidth', 3, 'Color', 'b', 'LineStyle', '-');   
    
    % negative slope
    yyaxis left
    id = maxidx(~sl(:, f), f);
    
    neg_subj = subj(~sl(:, f));
    neg_subj(isnan(id)) = [];

    neg_el = el(~sl(:, f));
    neg_el(isnan(id)) = [];

    id(isnan(id)) =[];
    neg_slope = betaInfo.x(f, id)./1000;
    h=histogram(neg_slope, fb{f}, 'FaceColor', 'r', 'FaceAlpha', 0.3, ...
        'EdgeColor', 'none', 'Normalization', 'probability'); hold on;
    
    yyaxis right
    [y_n, x_n] = ksdensity(neg_slope, 'BandWidth', b(f));
    plot(x_n, y_n, 'LineWidth', 3, 'Color', 'r', 'LineStyle', '-');   
        
    yyaxis left
    %xlim(xlims(f, :))
    yticks(0:10:40);
    xlabel('Max Beta Bin (kHz)');
    ylabel('Count');
    set(gca, 'FontSize', 13);
    yyaxis right
    set(gca,'YColor','none');
    box off;
    
    % mark boundary of formant range
    xline(betaInfo.x(f, 1)./1000 - h.BinWidth/2, ':k', 'LineWidth', 2)
    endidx = find(~isnan(betaInfo.x(f, :)), 1, 'last');
    xline(betaInfo.x(f, endidx)./1000 + h.BinWidth/2, ':k', 'LineWidth', 2)

    % using linear mized effects
    tbl = table();
    tbl.Maxb = [pos_slope, neg_slope]';
    tbl.Dir = [ones(length(pos_slope), 1); ones(length(neg_slope), 1)*2];
    tbl.subj = [pos_subj, neg_subj]';
    tbl.el = [pos_el; neg_el];
    
    lme = fitlme(tbl,'Maxb~Dir+(1|subj) + (1|el:subj)');
    disp(['-------------------F' num2str(f) '-------------------'])
    disp(lme);

    clear x* y*
end
    
% raster plot for showing beta 
figure('Renderer', 'painters', 'Position', [800 800 250*2 1*250]); 
subplot(1, 2, 1);
rast = cell2mat([inflections.b_f1]);
[~, k] = max(rast, [], 2);
[~, i] = sort(k);
imagesc(rast(i, :));
xlim([1 15]);
xticks([]); xlabel('Formant 1');
yticks([]); ylabel('Elelectrode');
caxis([0 0.03]);
colorbar;
set(gca, 'FontSize', 13);
title({'Sorted Beta Weights'});

subplot(1, 2, 2);
rast = cell2mat([inflections.b_f2]);
[~, k] = max(rast, [], 2);
[~, i] = sort(k);
imagesc(rast(i, :));
caxis([0 0.015]);
xticks([]); xlabel('Formant 2');
yticks([]); ylabel('Electrode');
set(gca, 'FontSize', 13);
colorbar;
bm = brewermap(18, 'BrBG');
bm = bm(5:18, :);
colormap(bm)

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections formant *model

%% Figure 2 (h, i) - F1 vs. F2 Slopes Over All Electrodes

formants = [1 2];
edgeColor = 'k';

figure('Renderer', 'painters', 'Position', [100 800 350*1 800]); 

subplot(2, 1, 1);
% find out max slope for normalization procedure
max_slope=[-Inf -Inf];
allel=0;
for s=1:length(SIDs)
    SID=SIDs{s};    
    if isfield(betaInfo, SID) && ~isempty(betaInfo.(SID).els) 

        % select all electrodes with significant linear or sigmoidal rsq
        elidx = sum(reshape(betaInfo.(SID).rsq([1 3], :, formants), ...
            [], length(betaInfo.(SID).els))>0, 'omitnan')>0;
        allel= allel+length(elidx);
        if sum(elidx, 'omitnan')>0
            max_slope(1)=max(max_slope(1), ...
                max(abs(betaInfo.(SID).beta_lin(elidx, formants(1), 1))));
            max_slope(2)=max(max_slope(2), ...
                max(abs(betaInfo.(SID).beta_lin(elidx, formants(2), 1))));
        end
    end
end

disp(['Number of electrodes = ' num2str(allel)])
all_slopes = [];

% scatter slope values
numel = 0; r=[];
for s = 1:length(SIDs)
    SID=SIDs{s};         
    if isfield(betaInfo, SID) && ~isempty(betaInfo.(SID).els)        
               
        % look for significant linear and sigmoidal fits
        rsq = squeeze(betaInfo.(SID).rsq([1 3], :, formants)); 
        
        if length(size(rsq))==3
            tmpidx=(rsq(2, :, 1)>0 | rsq(1, :, 1)>0)';                
            tmpidx(:, 2)=(rsq(2, :, 2)>0 | rsq(1, :, 2)>0)';
        else
            tmpidx=(rsq(2, 1)>0 | rsq(1, 1)>0)';                
            tmpidx(:, 2)=(rsq(2, 2)>0 | rsq(1, 2)>0)';
        end
        clear rsq        
                
        % get slope value from linear fit
        slope_pos=[betaInfo.(SID).beta_lin(:, formants(1), 1)>0 ...
            betaInfo.(SID).beta_lin(:, formants(2), 1)>0];  
        slope=[betaInfo.(SID).beta_lin(:, formants(1), 1) ...
            betaInfo.(SID).beta_lin(:, formants(2), 1)]; 
        
        idx = find(or(tmpidx(:, 1), tmpidx(:, 2)))';
                                               
        % mark electrodes with significant linear rsq in both formants
        color = tmpidx(idx, 1) & tmpidx(idx, 2);
           
        % z dim is either subject number or electrode number
        showSID = 1;
        if showSID
            z = repmat(str2double(SID(2:end)), 1, length(idx));
        else
            z = betaInfo.(SID).els(idx);
        end

        % slopes are normalized by maximum slope value
        sz = max([squeeze(betaInfo.(SID).rsq(1, idx, formants)) ...
            squeeze(betaInfo.(SID).rsq(3, idx, formants))], [], 2)*80;  % 0.8
        r=[r; sz./80];
        scatter3(slope(idx, 1)/abs(max_slope(1)), ...
            slope(idx, 2)/abs(max_slope(2)), z, sz, color, 'filled',...
            'HandleVisibility', 'off', 'MarkerFaceAlpha', 0.8, ...
            'MarkerEdgeColor', 'k', 'MarkerEdgeAlpha', 0.8); hold on;
                                                          
        clear lineColor f2 idx         
        
        numel = numel + sum(or(tmpidx(:, 1), tmpidx(:, 2)), 'omitnan');
        all_slopes = cat(1, all_slopes, slope./abs(max_slope));
    end     
end

xlabel(['F' num2str(formants(1)) ' slope']);
ylabel(['F' num2str(formants(2)) ' slope']);

% figure style
view(2); 
grid off
xticks(-1:0.5:1);
yticks(-1:0.5:1);
set(gca, 'FontSize', 15)

% formatting
yline(0, 'Color', 'k', 'HandleVisibility', 'off'); 
xline(0, 'Color', 'k', 'HandleVisibility', 'off');
h = refline(1, 0);
h.LineStyle='--';
h.HandleVisibility = 'off';
h.Color = 'k';
colormap(brewermap(2, 'PRGn'));
cbh = colorbar('southoutside');

% figure style
set(cbh,'YTick',[0.25 0.75], ...
    'TickLabels', {'one', 'both'});
clear cs cbh

% correlation 
[r, p] = corr(all_slopes(:, 1), all_slopes(:, 2));
disp(['Rho of slopes ' num2str(r) ' at p=' num2str(p)]);
disp(['Number of electrodes ' num2str(numel)]);

% linear mixed effect model
disp('----------------------- F1-F2 effects ----------------------- ')
tbl = table();
tbl.F1 = all_slopes(:, 1);
tbl.F2 = all_slopes(:, 2);
tbl.el = (1:length(all_slopes(:, 2)))';
tbl.subj = inflections.SID;
lme = fitlme(tbl,'F2~F1+(1|subj) + (1|el:subj)');
disp(lme);

% Figure 1 - F1 vs. F2 Direction of Slopes Table Over All Electrodes
grp = nan(2, height(inflections));
for f = 1:2
    var = ['_f' num2str(f)];
    grp(f, :) = (inflections.(['sl' var]) + 1) .* ...
        double(inflections.(['sr' var]) > 0 | ...
        inflections.(['lr' var]) > 0);
end
[tmptbl, ~, ~,LABELS] = crosstab(grp(1, :), grp(2, :));
tbl = tmptbl';
tbl([2 3], :) = tbl([3 2], :);
% ax = axes('Position',[0.55 0.65 .25 .25]);
ax = subplot(2, 1, 2);
imagesc(tbl, 'AlphaData', .2); hold on;
clear grp var

for x=1:size(tbl, 1)
    for y=1:size(tbl, 2)
        if num2str(tbl(x, y))>0
            text(y-0.15, x, num2str(tbl(x, y)), 'FontSize', 15);
        end
    end
end
xelements=str2double(LABELS(~cellfun(@isempty,LABELS(:, 2)), 2))+1;
yelements=str2double(LABELS(~cellfun(@isempty,LABELS(:, 1)), 1))+1;

% figure labels
labels_x={'n.s.', '-', '+'};
labels_y={'n.s.', '+', '-'};
xticklabels(labels_x(xelements));
set(gca, 'XTick', 1:5);
yticklabels(labels_y(yelements));
set(gca, 'YTick', 1:5);
xlabel(['F' num2str(formants(1)) 'slope sign']); 
ylabel(['F' num2str(formants(2)) 'slope sign']); 

% figure style
horzline(1.5, [], [], '-', 2.5);
vertline(1.5, [], [], '-', 2.5);
colormap(ax, brewermap(5, 'Greys'))

set(gca, 'FontSize', 15);

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections formant *model

%% Figure 2 - SUPP: F1 and F2 Tuning Direction on MNI Brain

formants = [1 2];

% initialize design electrode structure
for f = formants
    desel=struct();
    % 0 is negative, 1 is positive
    desel.conds = 0:1;
    desel.sz = repmat(55, 1, 2);
    desel.cols = brewermap(2, 'RdBu');
    desel.labels = {'negative', 'positive'};

    % extract F1/F2 slopes from inflections table
    conds = inflections.(['sl_f' num2str(f)]);  
    for s=1:length(SIDs)
        SID = SIDs{s};
        tblidx = inflections.SID==str2double(SID(2:end));
        desel.(SID).elid = inflections.el(tblidx);
        desel.(SID).condition=conds(tblidx);      
    end

    % colored by distance from the diagonal
    plotMNIElec(SIDs, desel, 'lh', [datapath '/pt_data']);
    title(['Formant ' num2str(f)]);
 
    % for figure formatting
    %print(fullfile([datapath 'figure_jpgs/mniBrain_LH_F' num2str(f) ...
    %         'dir.jpg']), '-djpeg', '-painters', '-r600');
    %cla;
    
    plotMNIElec(SIDs, desel, 'rh', [datapath '/pt_data']);
    title(['Formant ' num2str(f)]);
    
    % for figure formatting
    %print(fullfile([datapath 'figure_jpgs/mniBrain_RH_F' num2str(f) ...
    %         'dir.jpg']), '-djpeg', '-painters', '-r600');
    %cla;
end

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* inflections* ...
    formant *model 

%% Figure 2 - SUPP: Double peak, single peak

modelnames={'onset_maxDtL_maxDtLOnset_vowelOnset_F0_aud', ...
    'onset_maxDtL_maxDtLOnset_vowelOnset_F0_audF1F2'};

load('stim_info/mel_centerF.mat', 'binfrqs');


figure;
for s = SIDs
    SID = s{1};    
    corpusStrf=loadMultModelStrf(SID, modelnames, 'dimex', ...
        [datapath '/pt_data'], 1, '');
    if ~isempty(corpusStrf{1})
        rsq_aud = corpusStrf{1}.meanTestR.^2;
        
        for el = find(rsq_aud > 0.15) 
            tidx = find(inflections.SID == str2double(SID(2:end)) ...
                    & inflections.el == el);
            x = binfrqs(1:end-1) + diff(binfrqs)./2;
            if ~isempty(tidx)            
                featidx = 5:size(corpusStrf{1}.meanStrf, 1);
                [~, idx] = max(mean(corpusStrf{1}.meanStrf(featidx, ...
                    :, el)), [], 2);

                y = smooth(corpusStrf{1}.meanStrf(featidx, idx, el));
                if inflections.sl_f1(tidx) && ~inflections.sl_f2(tidx)
                    subplot(1, 2, 1);                   
                elseif ~inflections.sl_f1(tidx) && inflections.sl_f2(tidx)
                    subplot(1, 2, 2);                   
                end
                
                plot(x,y, 'k', 'LineWidth', rsq_aud(el)*4); hold on;
                ylabel('STRF Beta Weights (a.u.)');
                xlabel('Frequency (Hz, log scale)');
                yticks([]);
                xticks([]);
                set(gca, 'FontSize', 15)
                xlim([250 3000]);
                set(gca, 'XScale','log');
            end
        end
    end
end

subplot(1, 2, 1);
xline(800, '--k', 'LineWidth', 2);
title('Single Peak');

subplot(1, 2, 2);
xline(800, '--k', 'LineWidth', 2);
title('Double Peak');

clearvars -except *all subj *vow *SIDs datapath bef aft tps betaInfo* ...
    inflections formant *model

