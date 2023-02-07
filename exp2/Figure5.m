fig2 = figure;
exel = [3 54 212]; % example electrode numbers
%%  - - A/B - - stimulus schematic & effect schematic
ax(1) = subplot(1,2,1); cla; hold on;
si=1;
incltok =alldat(si).mean.f2>alldat(si).mean.f1;
gscatter(alldat(si).mean.f2(incltok),alldat(si).mean.f1(incltok), alldat(si).mean.inVowSpace(incltok),...
    [0.75 0.75 0.75;0 0 0]);%, 'filled');%, 'MarkerFaceColor', 'k');
plot(dimCont.M(1,2:end), dimCont.M(2,2:end),'-', 'color', 'k', 'LineWidth', 1);
rfl = refline(1,0); rfl.Color = [0.2 0.2 0.2]; rfl.LineStyle = ':'; rfl.LineWidth = 2;
set(gca, 'XDir', 'reverse'); set(gca, 'YDir', 'reverse')
xlabel('F2 (Hz)'); ylabel('F1 (Hz)');
% axis equal;
xlim([250 3000]); set(gca, 'XTickLabelRotation', 30);
ylim([180 1050]); yticks([250:250:1000]);
colormap(ax(1),rbcm);
legend('out of DIMEX range', 'in DIMEX range');
title('range of token formant values');

%% expected responses
si = 2;
% with zscore
f1 = alldat(si).mean.f1;
f2 = alldat(si).mean.f2;

protonames = {'zscore(f1(f1<f2))+zscore(f2(f1<f2))',  '-zscore(f1(f1<f2))-zscore(f2(f1<f2))',...
    'zscore(f1(f1<f2))-zscore(f2(f1<f2))', '-zscore(f1(f1<f2))+zscore(f2(f1<f2))',...
    'zscore(f1(f1<f2)).*zscore(f2(f1<f2))',  '-zscore(f1(f1<f2)).*zscore(f2(f1<f2))'};

protoTitles = {'F1 + F2', '-(F1 + F2)','F1 - F2', 'F2 - F1','F1*F2', '-F1*F2'};
cresp = cell(size(protonames));
for i = 1:length(protonames)
    cresp{i}  = nan(size(f1));
    cresp{i}(f1<f2) = eval(protonames{i});
    vowspMean(i,1:2) = grpstats(cresp{i}, alldat(si).mean.inVowSpace, 'nanmean');
    cresp{i} = reshape(cresp{i}, 10,10);
end
% subplot
spls = [3, 4, 7, 8, 11, 12];
for i= 1:length(protonames)
    protoax(i) = subplot(3,4, spls(i));
end

for i= 1:length(protonames)
    subplot(protoax(i));
    cla; hold on;
    imagesc('XData', alldat(si).f2u, 'YData', alldat(si).f1u, 'CData', cresp{i}');
    plot(dimCont.M(1,2:end), dimCont.M(2,2:end),'-', 'color', 'k', 'LineWidth', 1);
    colormap(protoax(i),pinkcurmap); xticks(alldat(si).f2u);yticks(alldat(si).f1u);
    ylabel('f1'); xlabel('f2');
    set(gca, 'XTickLabelRotation',60);
    set(gca, 'xdir', 'reverse');
    set(gca, 'ydir', 'reverse');
    axis tight;
    title(protoTitles{i});
    axis off;
end
cb1 = colorbar;
cb1.Location = 'east';
cb1.Position(1) = .49;
cb1.Label.String = 'response'; cb1.Label.Position(1)=-2;
cb1.Ticks = cb1.Limits; cb1.TickDirection ='out';
cb1.TickLabels = {'min', 'max'};

%% - - D - F single electrode panels

figure;
for i = 1:length(exel)
    exelax(i*3-1) = subplot(3,length(exel), 2*length(exel)+ i);
    exelax(i*3-2) = subplot(3,length(exel), i); % ero
    exelax(i*3) = subplot(3,length(exel),  length(exel)+i); % betas

end

% for each electrode: Responses by F1, Responses by F2, 2d color map of
% responses, bar plot of unique variance explained by F1, F2, IA effect
ccols = cool(10);
elsubj = linmAll.subj(exel);
maineffcol = autumn(3);
effcols = maineffcol([1 3 2],:);
effcols(4,:) = [0 0 0];

for cel = 1:length(exel)
    si = elsubj(cel);
    subjelnum = linmAll.el(exel(cel));
    xax = ((1:size(alldat(si).resp,2))-alldat(si).befaft(1))/100;

        %% main ERP only
    cresp = alldat(si).mean.resp(subjelnum,:,:);

    cax = subplot(exelax(cel*3-2));
    cla;hold on;
    f1gr = max(round((1:length(alldat(si).f1u))/3),1);
    for invi = 0:1
        ctr = alldat(si).mean.inVowSpace == invi;
        cpl(invi+1) = shadedErrorBar(xax*1000, nanmean(cresp(1, :, ctr),3), nansem(cresp(1, :, ctr),3), {'lineWidth', 2,'color', vowSpCol(invi+1,:)},1);
    end
    xlim([-100 600]); xticks([0 500]);
    vertline(0); horzline(0);
    ylim([-1 max(ylim)]);
    yticks([0 round(max(ylim))]);
    ylabel('HGA (z)');
    xlabel('time (s)');
    if cel == 1
        legend([cpl(1:2).mainLine],'out of vowel space', 'in vowel space');
    end
    title(sprintf('el %d, r^2:%0.2f, uv:%0.2f',...
        exel(cel), linmAll.rsq(exel(cel), 2),linmAll.rsq(exel(cel),7)- linmAll.rsq(exel(cel),2)  ))
   
    %% Effect betas and R^2 values
    plotuv = 2:4;
    npred = length(plotuv);
    cax = subplot(exelax(cel*3));
    cla;hold on;
    bm = 2;

    b=bar(1:3, linmAll.betas{bm}(exel(cel),plotuv),0.25,'FaceColor', 0.5*ones(1,3));
    b.FaceColor = 'flat';
    b.CData = effcols([1 2 4],:);

    text(1:3, -1*ones(size(plotuv)), num2str(round(linmAll.reduVar{2}(exel(cel),plotuv)'*100)), 'HorizontalAlignment', 'center')
    horzline(-0.9);

    ylim([-1.2 1]);yticks([-1:1]); yticklabels({'uR2', num2str((0:1)')});
    xlim([0.5 3.5]), xticks(1:length(plotuv)), xticklabels(linmAll.betaName{bm}(plotuv)); set(gca, 'XTickLabelRotation', 90);

    if cel == 1, ylabel('effect beta (a.u.)'); else,         ylabel('');     end

    %% two-d plot
    cax = subplot(exelax((cel-1)*3+2));
    cla;hold on;
    imagesc('XData', alldat(si).f2u, 'YData', alldat(si).f1u, 'CData', squeeze(alldat(si).mean.maxresp(subjelnum, :, :))');
    plot(dimCont.M(1,2:end), dimCont.M(2,2:end), 'k-');
    colormap(cax, pinkcurmap);
    cax.YDir = 'reverse';     cax.XDir = 'reverse';
    axis tight;
    xticks(alldat(si).f2u(1:2:end)); xticklabels(alldat(si).f2u(1:2:end));
    yticks(alldat(si).f1u(1:2:end)); yticklabels(alldat(si).f1u(1:2:end));
    xticks(xlim); xticklabels([min(alldat(si).f2u), max(alldat(si).f2u)]);
    yticks(ylim); yticklabels([min(alldat(si).f1u), max(alldat(si).f1u)]);
    xlabel('F2 (Hz)');ylabel('F1 (Hz)');
    cax.XTickLabelRotation=60;

end


%% print ex el figure
print(fullfile('figures','exeltask.eps'), '-depsc', '-painters')
%% load  single vowel dimex receptive fields
cfig = openfig('taskEl_dimexResp.fig');
for cel = 1:length(exel)
    cfig.Children(cel*2-1).Visible = 'off';
    cfig.Children(cel*2).XColor = 'none';
    cfig.Children(cel*2).YColor = 'none';
end
set(0, 'currentfigure', fig2);
for cel = 1:length(exel)
    subplot(5, 6, (cel+1)*6); cla; axis on;
    cax_copy = copyobj([cfig.Children((cel-1)*2+1:cel*2)],  fig2);
    subplot(5, 6, (cel+1)*6, cax_copy(2));
    colormap(gca, pinkcurmap);
end
for cel = 1:length(exel)
    subplot(5, 6, (cel+1)*6);
    cb=colorbar; cb.Visible = 'off';
end

