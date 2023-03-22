function [medoids, quantiles]=plotVariance(formants, vowels, cols, drawquant, ...
    plotfig, plotsilh)
% determine f1 and f2 ranges by sentence details 
        
    if nargin<5, plotfig=1; end
    if nargin<7, plotsilh=0; end 
    
    if plotfig && nargin<6
        figure;
    end
    
    % find unique vowels
    unique_vowels=unique(vowels);

    % plot scatter
    if nargin<3, cols = ametrine(length(english_vowels_all)); end
    if nargin<4, drawquant = 1; end
    
    allsilh = silhouette(formants(1:2, :)', vowels, 'cityblock');
    % permutation silhouette
    nperm = 500;
    permsilh = nan(length(vowels), nperm);
    
    file = 'formant_clusterPerm.mat';
    if exist(file, 'file')
        load(file);
    else        
        for p = 1:nperm
            disp(['Perm: ' num2str(p)])
            permsilh(:, p) = silhouette(formants(1:2, :)', ...
                vowels(randperm(length(vowels))), 'cityblock');
        end
        save(file, 'permsilh');
    end
%     disp(['Average Silhouette p-value: ' ...
%         num2str(1-((sum(mean(allsilh)>mean(permsilh)))/(nperm+1)))])
%     

    % color according to vowel transcription
    ctr=1;
    for j=1:length(unique_vowels)       
        isVow = cellfun(@(x)isequal(x,unique_vowels{j}),vowels);
        vowIdx = find(isVow);

        % prune formant space so not as dense, medoids make sense        
        if plotfig
            if ~isempty(vowIdx)
                scatter_col=repmat(cols(j, :), length(vowIdx), 1);
                ctr=ctr+1;
            end
        
            scatter(formants(2, vowIdx), formants(1, vowIdx), 25, scatter_col, 'filled', ...
            'HandleVisibility', 'off', 'MarkerFaceAlpha', 0.1); hold on;
            xlabel('F2 (Hz)'); ylabel('F1 (Hz)');
        end
        clear scatter_col         
    end
     
    % find medoids & mark
    quantile_range=[0.05 0.15 0.25 0.75 0.85 0.95];
    num_quantiles=length(quantile_range);
    quantiles=nan(length(unique_vowels),num_quantiles);
    medoids=nan(length(unique_vowels), 4);
    
    for j=1:length(unique_vowels)        
        vow_idx=cellfun(@(x)isequal(x,unique_vowels{j}),vowels);
        
        % so timit and dimex corpus color matches
        for f=1:4
            medoids(j, f)=quantile(formants(f,vow_idx), 0.5);
        end
        if plotfig
            scatter(medoids(j, 2), medoids(j, 1), 60, cols(j, :), 'filled', ...
                'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'k'); hold on;
            if plotsilh
                text(medoids(j, 2), medoids(j, 1), ...
                    num2str(mean(allsilh(vow_idx), 'omitnan'), 1), ...
                    'FontSize', 15, 'FontWeight', 'bold');
            end
        end
        
        if drawquant
            quantiles(j, :, 1) =quantile(formants(1,vow_idx),num_quantiles);
            quantiles(j, :, 2) =quantile(formants(2,vow_idx),num_quantiles);
            % plot quantile circlesx
            for q=1:num_quantiles/2        

                % calculate the radii from quantile position
                ra = mean([pdist([squeeze(medoids(j, 1:2));squeeze(quantiles(j, q, :))'], 'euclidean'), ...
                    pdist([squeeze(medoids(j, 1:2));squeeze(quantiles(j, q, :))'], 'euclidean')], 'omitnan');

                if plotfig
                    % draw ellipse               
                    ellipse(ra, ra, 0,medoids(j, 2),medoids(j, 1), cols(j, :)); hold on; 
                end
            end
        end

    end
    
    if plotfig
        xlim([450 3050]); ylim([100 1150]); % F1
        legend(unique_vowels);
        set(gca, 'xdir', 'reverse'); set(gca, 'ydir', 'reverse');
    else
        disp(['Median silhouette score over all vowel clusters: ' num2str(median(allsilh))]);
        
    end
    
end