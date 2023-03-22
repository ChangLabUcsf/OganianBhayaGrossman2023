% perform decoding of labels and permutes labels, reps time
function [acc, permacc, beta, wind, nanel, trls, classes] ...
    = decode_permute(X, y_fstat, y, trialExcl, reps, pcaflag)
    % X - matrix of form electrodes by timepoints by trials
    % y - category labels
    % trialExcl - a scalar (0, 1) to determine how often to exclude trials, 
    % higher values remove fewer trials and more electrodes, lower values
    % exclude more trials
    % reps - number of repetitions
    
    % for Fstat, remove all nanidx trials
    % remove columns and rows where data is missing
    % trials for which over electrodes do not have data
    % remaining electrodes with missing trials
    n = 100;
    if isempty(trialExcl)
        % find trialExcl param that retains at least 100 trials per class
        trialExcl = linspace(0, 1, n);
        
        trls = ~(squeeze(sum(isnan(mean(X, 2)), 1) > trialExcl*size(X, 1)));
        mintrls = nan(n, 1);
        minel = nan(n, 1);
        for i = 1:n
            mintrls(i) = min(crosstab(y_fstat(trls(i, :)))); 
            minel(i) = sum(~(sum(isnan(squeeze(mean(X(:, :, trls(i, :)), 2))), 2)>0));
        end
        
        debugFig = 0;
        if debugFig
            figure;
            scatter(mintrls, minel,55, trialExcl, 'filled');
            ylabel('Minimum number of trials per category');
            ylabel('Electrodes');
        end
        trialExcl = trialExcl(find(mintrls>100, 1));
    end
    
    nantrls = squeeze(sum(isnan(mean(X, 2)), 1) > trialExcl*size(X, 1));
    X(:, :, nantrls) = [];
    nanel = sum(isnan(squeeze(mean(X, 2))), 2)>0;
    X(nanel, :, :) = [];
    
    % median time point for Fstat across electrodes
    numclass = length(unique(y_fstat(~nantrls)));
    [Fstat, ~, ~, ~, ~] = Fstat_TIMIT(X, y_fstat(~nantrls), ...
        unique(y_fstat(~nantrls)));          
    [~, maxtp]=max(Fstat, [], 2);    
    wind = round(median(maxtp)-2:median(maxtp)+2);
    X = X(:, wind, :);

    acc = cell(reps, 1);
    permacc = cell(reps, 1);
    
    % minimum electrodes
    numel = size(X, 1);
    warning(['Final electrode set: ' num2str(numel)]); 
    warning(['Min number of trials in category: ' num2str(min(crosstab(y(~nantrls))))])
    warning(['Window used:' num2str(wind)]);
       
    % find which trials retained
    trls = find(~nantrls);

    % run SVM decoder on all time points
    data = struct();
    data.resp = X;
    data.vowel = y(~nantrls);

    % beta is of form pair x electrode x rep
    % always uses a binary classifier
    beta = nan(nchoosek(numclass, 2), numel, reps);
    for rep = 1:reps
        
        % nonpermuted data case
        [svmOut, tmp,~, ~, ~] = ecog_svm_main(data, 1:size(X, 1), ...
                                    'resp', 'vowel', 1:length(trls), ...
                                    '', 1:size(X, 2), pcaflag, 1);

        % save out beta weights per electrode (or PC) per classifier
        weights = cell2mat({svmOut.beta}');
        beta(:, 1:size(weights, 2), rep) = weights;
        
        a = tmp;
        a(isnan(tmp)) = [];
        acc{rep} = a;

        % permuted data case
        data_perm = data;
        idx = randperm(length(trls));
        data_perm.vowel = data_perm.vowel(idx);
        [~, tmp,~, ~, ~] = ecog_svm_main(data_perm, 1:size(X, 1), ...
                                    'resp', 'vowel', 1:size(data.resp, 3), ...
                                    '', 1:size(X, 2), pcaflag, 1);
        a = tmp;
        a(isnan(tmp)) = [];
        permacc{rep} = a;
        clear a data_perm
    end

    classes = {svmOut.classes};
end