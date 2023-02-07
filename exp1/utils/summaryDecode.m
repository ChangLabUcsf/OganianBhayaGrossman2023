function summaryDecode(Dvow, acc_all, pthresh)
    
    load('svmClassOrder.mat', 'classOrd');
    ord = [classOrd.classesNum];

    % plot all the medoids
    medoids = plotVariance(Dvow.formantVals, Dvow.vowel, [], 0, 0);
    
    figure;
    scatter(medoids(:, 2), medoids(:, 1), 95, ...
        getColors(1), 'filled', 'MarkerEdgeColor', 'k'); 
    hold on;
    
    % line from medoid to medoid (color is difference, thickness is magnitude)
    % dashed is non-significant?
    ctr = 1;
    cmax = max(abs(squeeze(diff(median(acc_all(:, 1:2, :), 1)))));
    for pair = ord
        v1 = pair(1);
        v2 = pair(2);
        acc = squeeze(acc_all(:, :, ctr));
            
        % determine significance vs permutation distribution
        % h, p(i)] = ttest2(decode_all.perm_acc(pair, :, i), acc(:, i)); 
        [~, p] = ttest2(acc(:, 1), acc(:, 2)); 
        disp(num2str(p));

        % determine magnitude
        n = 7;
               
        cspec = brewermap(n, 'RdBu');
        cspec(ceil(n/2), :) = [0.7 0.7 0.7];

        cidx = discretize(diff(median(acc(:, 1:2))), linspace(-1*cmax, cmax, n));
        width = discretize(max(median(acc(:, 1:2))), linspace(0.45, 0.75, 5))*0.75;

        style = '-';
        if p>pthresh, style = ':'; end

        plot([medoids(v1+1, 2) medoids(v2+1, 2)], ...
            [medoids(v1+1, 1) medoids(v2+1, 1)], 'Color', cspec(cidx, :), ...
            'LineWidth', width, 'LineStyle', style);
        ctr = ctr+1;
    end
    
    % formatting
    set(gca, 'YDir', 'Reverse', 'XDir', 'Reverse');
    ylim([300 700]);
    xlim([1000 2500]);
    yticks(300:100:800);
    xticks(1000:500:2500);
    yticklabels(split(num2str((300:100:700)./1000)));
    xticklabels(split(num2str((1000:500:2500)./1000)));
        
    set(gca, 'FontSize', 15);
    ylabel('F1 (kHz)');
    xlabel('F2 (kHz)');
    
    % for colorbar
    figure;
    scatter(1:n, 1:n, 95, ...
        cspec, 'filled');
    colormap(cspec);
    caxis([-0.15, 0.15]);
    cbh = colorbar;
    set(cbh, 'ticks', [-0.15, 0.15]);
    set(gca, 'FontSize', 15);
    ylabel(cbh, 'diff(F1-/F2+, F1+/F2-)');   
end