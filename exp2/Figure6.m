%% ----- -----  Figure 6: Population analysis ----- ----- -----


%% A - - boxplot of unique % variance explained by each predictor
figure,
subplot(2,2, [1 3]); cla; hold on;
plotuv = 2:6;
boxplot(linmAll.reduVar{2}(AnaElTask, 2:4), ...
    'BoxStyle', 'filled', 'PlotStyle', 'compact', 'Labels', linmAll.betaName{2}(2:4), 'Colors', 'k');
yticks(0:.2:.8);
ylabel('unique R2');horzline(0);
title('Unique R2 by effect');
%% - -  - - Variance partitioning in model including full stimulus set
plotuv = [2 3;1 4]; % 3 plots showing different combinations of unique variances
elsizeEH = max(max(linmAll.reduVar{2}(:,2:3)*100,[],2),10);
IAEl = linmAll.reduVar{2}(:,4) > 0.05;
UVvals = [linmAll.rsq(:,1), linmAll.reduVar{1}(:,2:3), linmAll.reduVar{2}(:,4), max(linmAll.reduVar{1}(:,2:3),[],2)];
cuvnames = { 'F1+F2','F1', 'F2', 'F1*F2', 'max(F1,F2)'};
%% B+C - - Unique variance scatter plots

for i = 1:size(plotuv,1)
    subplot(2,2,i*2); cla;hold on;
    cx = UVvals(AnaElTask,plotuv(i,1));
    cy=UVvals(AnaElTask,plotuv(i,2));
    [r,p]=corr(cx,cy)
    scatter(cx, cy, 15,'k', 'filled');
    xlabel(['R2(' cuvnames{plotuv(i,1)} ')']), ylabel(['R2(' cuvnames{plotuv(i,2)} ')']);
    rfl = refline(1,0); rfl.Color = 'k'; rfl.LineStyle =':';
    [cr,cp]=corr(cx, cy);
    ylim([-0.1 max(cy)*1.1]);
    text(0.1, max(ylim)*.95, sprintf('r = %0.2f, p = %2.0e', cr,cp));
    horzline(0,[],'k', '-',1); vertline(0,[],'k', '-',1);
    xticks([0:.4:.8]); yticks([0:.2:.6]);
end

%% unique variance correlation with mixed effects model
% model: uv2 ~ uv1 + (1 | el:SID)
clme = cell(2,1);
for i = 1:size(plotuv,1)
    cx = UVvals(AnaElTask,plotuv(i,1));
    cy=UVvals(AnaElTask,plotuv(i,2));
    csubj= linmAll.subj(AnaElTask)';
    [r,p]=corr(cx,cy)
    lmetable = table(cx, cy, csubj); 
    clme{i} = fitlme(lmetable, 'cy ~ cx + (1|csubj)+(1|csubj:cx)');
end
%% - - - D - - - betas 
clear cpl;
figure;
subplot(2,2,[1 2]); cla; hold on;

maineffcol = autumn(3);
effcols = maineffcol([1 3 2],:);
effcols(4,:) = [0 0 0];
effsizes = [70, 70, 70,30,30];
effmarker = 'ooo^v';


rsqCutOff=0;
plotel=cell(4,1);
% only F1
plotel{1} = (UVvals(:,2) >rsqCutOff & UVvals(:,3) < rsqCutOff &  AnaElTask);
% only F2
plotel{2} = (UVvals(:,3) > rsqCutOff & UVvals(:,2) < rsqCutOff & AnaElTask);
% F1 and F2
plotel{3} = (UVvals(:,2) > rsqCutOff & UVvals(:,3) > rsqCutOff & AnaElTask);
% linear Interaction
plotel{4} = (UVvals(:,4) > rsqCutOff & AnaElTask);

% plot
for cplid = 1:3
    cpl(cplid)=scatter(linmAll.betas{2}(plotel{cplid}, 3), linmAll.betas{2}(plotel{cplid}, 2),...
        effsizes(cplid),effcols(cplid,:),'filled', 'Marker', effmarker(cplid), 'MarkerEdgeColor','k'); % rcm(elcols(plotel),:), 'filled');
end
cplid = 4;
cpl(cplid)=scatter(linmAll.betas{2}(plotel{cplid}, 3), linmAll.betas{2}(plotel{cplid}, 2),...
    effsizes(cplid),effcols(cplid,:),'filled', 'Marker', effmarker(cplid)); % rcm(elcols(plotel),:), 'filled');

rfl(1) = refline(1,0);
rfl(2) =refline(-1,0);
for i = 1:2, rfl(i).Color = 'k'; rfl(i).LineStyle =':';end

% axis tight;
xlim([-1 1]); ylim([-1 1]);
horzline(0); vertline(0);
xlabel('F2 beta'); ylabel('F1 beta')
set(gca, 'XDir', 'reverse');set(gca, 'YDir', 'reverse');
title('Main Effect beta values');
legend(cpl,{'F1', 'F2', 'F1&F2', 'F1*F2'});


% correlation
[cr,cp]=corr(linmAll.betas{2}(AnaElTask, 3), linmAll.betas{2}(AnaElTask, 2));
% mixed effects model
cx = linmAll.betas{2}(AnaElTask, 3);
cy=linmAll.betas{2}(AnaElTask, 2);
csubj= linmAll.subj(AnaElTask)';
lmetable = table(cx, cy, csubj);
clme = fitlme(lmetable, 'cy ~ cx + (1|csubj)+(1|csubj:cx)');

%% count of electrodes by  tuning direction
crosstab(sign(linmAll.betas{2}(AnaElTask, 3)), sign(linmAll.betas{2}(AnaElTask, 2)))

%% - - - D 

% add to plot
text(0.1, min(ylim)*.95, sprintf('r = %0.2f, p = %2.0e', cr,cp));

%% comparison to model in vowel space
cmod=2;
subplot(2,2,3)
cla;hold on;
i=1;
cols = 'rb';
for cs1 = [-1 1]
    for cs2 = [-1 1]
        cel = AnaElTask & sign(linmAll.betas{2}(:,2))==cs1  & sign(linmAll.betas{2}(:,3))==cs2;
        scatter(linmAll.rsq(cel,2),linmAllSm.rsq(cel,2), 20, cols((cs1==cs2)+1), 'filled');
        i=i+1;
    end
end
xlim([0 1]); ylim([-1 1]);
yticks(-1:.5:1);
refline(1,0); horzline(0); vertline(0);
legend('on-diagonal', 'off-diagonal')
xlabel('R^2 full model'), ylabel('R^2 vowel space model');
xticks([0 .5 1]); yticks([-1:.5:1]);
subplot(2,2,4)
cla; hold on;
rdiff = linmAll.rsq(AnaElTask,2)-linmAllSm.rsq(AnaElTask,2);
boxplot(rdiff, [sign(linmAll.betas{2}(AnaElTask,2)), sign(linmAll.betas{2}(AnaElTask,3))],...
    'Labels', {'-1/-1', '-1/1', '1/-1', '1/1'},'Positions', [1 3 4 2], 'PlotStyle', 'compact', 'Colors', 'k');
horzline(mean(linmAll.rsq(AnaElTask,cmod)-linmAllSm.rsq(AnaElTask,cmod)));
title('R^2(full space) - R^2(vowel space)');
yticks(-1:.5:1);

%% stats

offdiag = sign(linmAll.betas{cmod}(:,2))== sign(linmAll.betas{cmod}(:,3));
[~,p1,~,stat1]=ttest(linmAll.rsq(AnaElTask,cmod)-linmAllSm.rsq(AnaElTask,cmod));

[~,p1,~,stat1]=ttest(linmAll.rsq(offdiag & AnaElTask,cmod)-linmAllSm.rsq(offdiag & AnaElTask,cmod));
[~,p2,~,stat2]=ttest(linmAll.rsq(offdiag==0 & AnaElTask,cmod)-linmAllSm.rsq(offdiag==0 & AnaElTask,cmod));

%% mixed effects stats

lmemat = [offdiag(AnaElTask), linmAll.subj(AnaElTask)', linmAll.rsq(AnaElTask,cmod), (1:sum(AnaElTask))', ones(sum(AnaElTask),1)];
lmemat2 = [offdiag(AnaElTask), linmAll.subj(AnaElTask)', linmAllSm.rsq(AnaElTask,cmod), (1:sum(AnaElTask))', 2*ones(sum(AnaElTask),1)];
lmetable = array2table([lmemat;lmemat2]);
lmetable.Properties.VariableNames= {'offdiag', 'subj', 'cx', 'el', 'mod'} ;

fitlme(lmetable, 'cx ~mod*offdiag+(1|subj) + (1|el:mod)')


