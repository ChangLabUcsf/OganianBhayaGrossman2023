function [useFigure] = ecog_erpPlotSingle(out, desMat,useEl,gridLayout ,plotTimePoints,useFigure, bef, onsettimes, plotMean, dataf,elnames, colors)

% function [useFigure] = ecog_erpPlot(out, desMat,useEl,gridLayout ,plotTimePoints,useFigure, bef, onsettimes, plotMean, dataf)
% plot averag ERPs for data in out structure by conditions indicated in conditions
% out - out structure or data matrix of form el x time x trial
% desMat: structure with two fields. desMat.condition is a vector of the same
% length as out structure/matrix(3rd D) with numbers from 1 to # different plot
% conditions and desMat.names contains condition names for each of them
% returns handle to figure
% bef - how much time (in sec) before stim. onset?
% onsettimes - tp to which to align every single trial
% plotMean 0 plot mean erp (1) or single trials (0)
% (c) Yulia Oganian, July 2016

%% argin
if nargin < 9, plotMean = 1; end
if nargin < 10, dataf= 100; end
if nargin < 5, plotTimePoints = []; end
if nargin < 11 elnames=[];end

%% colors for plot
if nargin < 12
    colors = ametrine(3); % morgenstemning , ametrine
    if length(desMat.names)>3
        colors = ametrine(length(desMat.names));
    end
    %colors = flip(colors);
    %colors = colors(end-5:end, :);
    %colors = colors(3:3:27, :);
end

if nargin < 6 || isempty(useFigure)% use new figure
    new_fig = 1;
    useFigure  =figure;
    nCondSoFar = 0;
else
    %figure(useFigure); % plot to existing figure
    new_fig = 0;
    %% choose colors for additional condition
    l = legend('Show');
    condNamesSoFar = l.String;
    nCondSoFar = length(condNamesSoFar);
    usedColors = nan(nCondSoFar,3);
    for i = 1:nCondSoFar
        chandles = findobj(gcf, 'DisplayName', l.String{i});
        lineHandles(i)  = chandles(1);
        usedColors(i,:) = lineHandles(i).Color;
    end
    
%     colors = colors(~ismember(colors,usedColors, 'rows'),:);
%     colors = colors(end:-1:1,:);
    hold on;
end
if isempty(useEl), useEl = gridLayout(:);end
if nargin < 7, bef = .5; end
if isempty(bef), bef = .5; end
if nargin < 8, onsettimes = []; end

%% prepare data
if isempty(onsettimes) % all trials already time-aligned
    if isstruct(out) % out structure
        a = {out.resp}';
        sizes = cellfun('length', a);
        longestSent = max(sizes);
        allresp = nan(size(a{1},1),longestSent, length(a));
        for i = 1:length(a) % mean across presentation of the same item
            allresp(:,1:size(a{i},2),i) = mean(a{i},3);
        end
        %         dataf = out(1).dataf;
    elseif isnumeric(out)
        allresp = out;
        %         dataf = 100; % assume standard data frequency of 100 Hz
    end
else % need to realign
    if isstruct(out) % out structure
        disp('no realignement procedure for out structure.')
        return;
    elseif isnumeric(out)
        dataf = 100; % assume standard data frequency of 100 Hz
        newrespDur = size(out,2)-ceil(max(onsettimes)*dataf);
        allresp = nan(size(out,1),newrespDur,size(out,3));
        for i = 1:size(allresp,3)
            allresp(:,:,i) = out(:, floor(onsettimes(i)*dataf):(newrespDur+floor(onsettimes(i)*dataf)-1) ,i);
        end
    end
end

if isempty(plotTimePoints), plotTimePoints = 1:size(allresp,2); end

if isempty(desMat)
    desMat.condition = ones(size(allresp, 3),1);
    desMat.names = {'mean HGA'};
else
    if ~isfield(desMat, 'condition') || isempty(desMat.condition), desMat.condition = ones(size(allresp, 3),1);end
    if ~isfield(desMat, 'names')|| isempty(desMat.names), desMat.names = {'mean HGA'}; end
end
nNewCond = length(desMat.names);
%% calculate mean and sem for each condition
nCond = length(desMat.names);
condNums = unique(desMat.condition);
if length(condNums) ~=nCond
    fprintf(2,'wrong number of conditions. check design matrix. exiting.');
    return
end
meanErps = nan(size(allresp,1), length(plotTimePoints), nCond);
semErps = nan(size(meanErps));
for ccond = 1:length(desMat.names)
    meanErps(:,:,ccond) = mean(allresp(:,plotTimePoints,desMat.condition == condNums(ccond)), 3, 'omitnan');
    semErps(:,:,ccond) = nansem(allresp(:,plotTimePoints,desMat.condition == condNums(ccond)), 3);
end
%% choose subgrid to plot

[row, col] = find(ismember(gridLayout, useEl));
curGrid = gridLayout(min(row):max(row), min(col):max(col));
%% plot mean yes/no
if plotMean == 1
    temperps = meanErps(useEl,:,:);
    ymin = min(temperps(:))*1.2;
    ymax = max(temperps(:))*1.2;
else
    ymin = min(reshape(allresp(useEl,plotTimePoints,:),[],1,1))*1.2;
    ymax = max(reshape(allresp(useEl,plotTimePoints,:),[],1,1))*1.2;
end

time_axis= (0:(size(allresp,2)-1))/dataf - bef;
time_axis = time_axis(plotTimePoints);


%% plot grand average for each electode separately
elcount = 1;
plotSEM = 1;
lineWidth = 2;
for i = 1:numel(curGrid)
    cel = curGrid(i);
    fprintf(2, 'plotting el %d out of %d ... \n', i, numel(curGrid))
    if ~isnan(cel)
        %subplot('Position',p);
        hold on
        for ccond = 1:nCond
            if plotMean
                if plotSEM
                    chl(ccond) = shadedErrorBar(time_axis,...
                        meanErps(cel, :, ccond),...
                        semErps(cel, :, ccond), 'lineprops', {'color', colors(ccond,:), 'LineWidth', lineWidth});
                else
                    chl(ccond) = shadedErrorBar(time_axis,...
                        meanErps(cel, :, ccond),...
                        zeros(size(meanErps(cel, :, ccond))), 'lineprops', {'Color', colors(ccond,:), 'LineWidth', lineWidth});
                end
            else
                for cpl = 1:size(allresp,3)
                    plot(time_axis, squeeze(allresp(cel,:,cpl)), 'color', colors(ccond,:), 'LineWidth', lineWidth);
                end
            end
        end
        
        %% background color
        % color background of plot if electrode was supposed to be plotted
        % & electrode is not the only one being plotted
        if ~ismember(cel, useEl) && size(useEl, 1) ~= 1
            set(gca,'Color',[0.8 0.8 0.8]);
        end
        axis tight;
        line(get(gca,'XLim'),[0 0],'Color','k');
        if i ~= 1
            set(gca,'XTickLabel',[]);
            set(gca,'YTickLabel',[]);
        end
        if i == 1 %adjust ylim
            if ~new_fig
                ylimold = get(gca, 'YLim');
                ymin = min(ylimold(1), ymin);
                ymax = max(ylimold(2), ymax);
            end
        end
        ylim([ymin,ymax])
        line(get(gca,'XLim'),[0 0],'Color','k');
        line([0 0], get(gca,'YLim'),'Color','k');
        
        % el numbers
        text(min(get(gca, 'XLim'))*.8,max(get(gca,'YLim'))*0.8,num2str(cel), 'Color', 'r');
        if mod(i, size(curGrid, 1))==1
            if ~isempty(elnames)
                title(elnames{elcount}(1:end-1))
            end
        end
        elcount = elcount +1;

        %% legend
        if i == 1
            for ccond = 1:nCond
                lineHandles(nCondSoFar + ccond) =  plot(time_axis,...
                    meanErps(cel, :, ccond), 'color', colors(ccond,:));
            end
            % legend
            if new_fig == 1
                legend([chl.mainLine],desMat.names);
            else
                legend([lineHandles [chl.mainLine]], [condNamesSoFar, desMat.names])
            end
        end
    end
    alpha(0.3)
end

