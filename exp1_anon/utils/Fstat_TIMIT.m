function [fvals, betweenVar, withinVar, df1, df2] = Fstat_TIMIT(D, L, phns)
% function [Fstat, betweenVar, withinVar, df1, df2] = Fstat_TIMIT(D, L, phns)
%   
% Input:
%   D: matrix with [electrodes x time x trials]
%   L: vector of length [trials], each number is the phoneme class (int)
%   phns: vector of numbers representing the unique classes
%   [unique(trials)]  (ex: phns = [1:33])
%
% Output:
%   Fstat: f statistic [electrodes x time] 
%   betweenVar: Between phoneme variance [electrodes x time]
%   withinVar: within phoneme variance [electrodes x time]
%   df1: degrees of freedom for numerator
%   df2: degrees of freedom for denominator
%
% Liberty Hamilton 2016
%

overallMean = nanmean(D,3); % mean across all phonemes

% Initialize variables
withinSS = zeros(size(D,1),size(D,2));
betweenSS = zeros(size(D,1),size(D,2));
fvals = zeros(size(D,1),size(D,2));

nphn = length(phns); % Get total number of classes
total_trials = 0; % initialize trial #

% Loop over all phoneme classes
for phoneme = 1:nphn
    
    Dphn = D(:,:,L==phns(phoneme)); % find all trials with this particular phoneme
    
    withinMean = mean(Dphn,3); % mean across all trials for this phoneme
    
    % Within class sum of squares
%     withinSS = withinSS + sum((Dphn-repmat(withinMean,[1,1,size(Dphn,3)])).^2,3);
    withinSS = withinSS + sum((Dphn-withinMean).^2,3);
    % Between class sum of squares
    betweenSS = betweenSS + size(Dphn,3).*(mean(Dphn,3) - overallMean).^2;
    
    % Increment trial #
    total_trials = total_trials + size(Dphn,3);
end

% Calculate the degrees of freedom
df1 = nphn - 1; % degrees of freedom for numerator
df2 = total_trials - nphn; % degrees of freedom for denominator

% Get between phoneme and within phoneme variance estimates
betweenVar = betweenSS/df1; 
withinVar = withinSS/df2;

fvals = betweenVar./withinVar;
