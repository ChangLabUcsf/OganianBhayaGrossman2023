
addpath(genpath('utils'));
datapath = '../data/';

% read in relevant CSV columns
file = [datapath 'mturk_results/Batch_3.csv'];      
T = readtable(file, 'Delimiter', ',');
vowel_types = {'[a] as in had', '[i] as in hid', '[ee] as in heed', ...
     '[o] as in hop', '[u] as in hug', '[oo] as in hoop', '', 'Does not sound like a vowel'};
T(strcmp(T.Answer_word, 'Cannot determine'), :) = [];

% read in F1-F2 location, whether perceived as human, which vowel
url_split = cellfun(@(x) split(x, '_'), T.Input_audio_url, ...
    'UniformOutput', false);
T.F1 = cellfun(@(x) str2double(x{4}), url_split);
T.F2 = cellfun(@(x) str2double(x{6}(1:end-4)), url_split);
T.human = cellfun(@(x) ismember({'Yes'}, x), T.Answer_word);
T.vowel = cellfun(@(x) find(ismember(vowel_types, x)), T.Answer_vowel);

% plot the percentage of responses that said synthesized token sounds like
% "human"/which vowel
human_plot = [];
vowel_plot = [];
for f1 = unique(T.F1)'
    for f2 = unique(T.F2)'
        idx = T.F1 == f1 & T.F2 == f2;
        if ~isempty(idx)
            human_plot = [human_plot; f1, f2, mean(double(T.human(idx))), ...
                sum(ismember(T.vowel(idx), 8))];
            vowel_plot = [vowel_plot; f1, f2, sum(ismember(T.vowel(idx), [1, 5])), ... % combine /a, uh/
                sum(ismember(T.vowel(idx), [2, 3])), ... % combine /i, e/
                sum(ismember(T.vowel(idx), 6))]; % combine /o, oo/
        end
    end
end

% plot showing the percentage of responses saying "human produced"
figure;
subplot(1, 2, 1);
scatter(human_plot(:, 2), human_plot(:, 1), 35, ...
    human_plot(:, 3), 'filled'); hold on;

stimdir = '../exp2/data/stim';
dimCont = load(fullfile(stimdir,'dimexContour.mat'));
[dimCont.M, c] = contour(dimCont.C{1}, dimCont.C{2},dimCont.N', [50, 50]);
c.LineColor='k';
colormap(brewermap(200, 'RdBu'));
caxis([0.5 1]);
colorbar()
vowelFormat()

% comparing "human produced" percentage outside and inside formant space
subplot(1, 2, 2);
in = inpolygon(human_plot(:, 2),human_plot(:, 1),...
    dimCont.M(1, 2:end), dimCont.M(2, 2:end));
boxchart(ones(sum(in), 1), human_plot(in, 3), ...
    'BoxFaceColor',[0.5 0.5 0.5], 'MarkerStyle','+', 'MarkerColor','k'); hold on;
boxchart(2*ones(sum(~in), 1), human_plot(~in, 3), ...
    'BoxFaceColor',[0.5 0.5 0.5], 'MarkerStyle','+', 'MarkerColor','k');
[h, p, ci, stats] = ttest2(human_plot(in, 3), human_plot(~in, 3));
title(['P-Value: ' num2str(p)]);

ylim([0.3, 1.1]);
yticks(0.5:0.5:1)
xticks([1, 2]);
xticklabels({'Inside', 'Outside'});
xlabel('Relative to natural formant space');
set(gca, 'FontSize', 13);

% looking at vowel identity responses across tokens
figure;
subplot(1, 3, 1)
[cols] = getColors(1);
cols = cols([1, 3, 5], :);
for i = 1:3
    ax = subplot(1, 3, i);
    [dimCont.M, c] = contour(dimCont.C{1}, dimCont.C{2},dimCont.N', [50, 50]); hold on;
    c.LineColor='k';

    scatter(vowel_plot(:, 2), vowel_plot(:, 1), 35, vowel_plot(:, 2+i), 'filled');
    colormap(ax, [linspace(1, cols(i, 1)); ...
        linspace(1, cols(i, 2)); linspace(1, cols(i, 3))]');
    colorbar
    caxis([0 9]);
    vowelFormat()
end


function vowelFormat()
    set(gca, 'Xdir', 'reverse', 'YDir', 'reverse');    

    xticks([500 3000]);
    xticklabels({'0.5', '3'});   
    xlabel('F2 (kHz)');
    xlim([500 3000]);

    yticks([200, 1000]);
    yticklabels({'0.2', '1'});
    ylabel('F1 (kHz)');
    ylim([200, 1000]);
   
    set(gca, 'FontSize', 13);
end

