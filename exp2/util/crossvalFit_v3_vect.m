%% cross-validated fitted curve _v3
function [y_pred, param, rsq, RSS, rCorr,Ftest,redMod]=crossvalFit_v3_vect(x, y, type, numBins)

if nargin<4, numBins = length(x); end

nchan = size(y,2);
% cross validation is leave one out
y_pred = nan(numBins,nchan);

% coefficients of cross-validated fitted curves
switch type
    case 'sigmoid'
        param = nan(3  + size(x,2), nchan,numBins);
    case 'linear'
        param = nan(nchan,size(x,2),numBins);
        param_red = cell(size(x,2),1);
        y_pred_red = cell(size(x,2),1);
        for i = 1:size(x,2)
            param_red{i} = zeros(size(param));
            y_pred_red{i} = zeros(size(y_pred));
        end
end
%% leave-one-out cv
for excl = 1:numBins
    % use subset of x, and y to fit curve
    idx = (1:size(x,1)) ~=excl;
    x_excl = x(idx,:);
    y_excl = y(idx,:);
    pidx = 1:size(x,2);
    switch type
        case 'linear'
            % full model
            param(:,:,excl) = glm_vect(x_excl, y_excl);
            y_pred(excl,:) = param(:,:,excl)*x(excl,:)';
            
            % reduced models, leave one parameter out
            for i = 1:size(x_excl, 2)
                param_red{i}(:,pidx~=i,excl) = glm_vect(x_excl(:, pidx~=i), y_excl);
                y_pred_red{i}(excl,:) = param_red{i}(:,:,excl)*x(excl,:)';
            end
        case 'sigmoid'
            for cchan  = 1:size(y,2)
                
                
                param(:,cchan, excl) = sigm_fit_multx(x_excl, y_excl(:,cchan), [], ...
                    [min(y_excl(:,cchan)) max(y_excl(:,cchan)) 0 zeros(1,size(x,2))], 0);
                fsigm = @(param,xval) ...
                    param(1)+(param(2)-param(1))./(1+10.^(param(3)-param(4:end)*xval'));
                y_pred(excl,cchan) = fsigm(param(:,cchan, excl)', x(excl,:));
                
                % reduced models, leave one parameter out
                for i = 1:size(x_excl, 2)
                    param_red{i}([1 2 3 find(pidx~=i)+3],cchan,excl) = sigm_fit_multx(x_excl(:,pidx~=i), y_excl(:,cchan), [], ...
                        [min(y_excl) max(y_excl) 0 zeros(1,size(x,2)-1)], 0);
                    fsigm = @(param,xval) ...
                        param(1)+(param(2)-param(1))./(1+10.^(param(3)-param(4:end)*xval'));
                    y_pred_red{i}(excl,cchan) = fsigm(param(:,cchan, excl)', x(excl,:));
                end
            end
            param = permute(param, [2 1 3]);
    end
end

% full model rsq value can be negative based on how this is calculated
RSS = sum((y - y_pred).^2,1);
SStotal = sum((y - mean(y,1)).^2);
rsq = 1-(RSS./SStotal);
rCorr = (diag(corr(y,y_pred)))';


% reduced models
for i = 1:size(x_excl,2)
    redMod.RSS_red(i,:) = sum((y - y_pred_red{i}).^2);
end
redMod.rsq_red = 1-(redMod.RSS_red./SStotal);

%% unique variance
redMod.uvar = rsq - redMod.rsq_red;

%% comparison Ftest
fdiff = @(rss1,rss2,p1,p2,n) ((rss1-rss2)./(p2-p1))./(rss2./(n-p2));
Ftest.df = [1, length(y)- size(x,2)];

for i = 1:size(redMod.RSS_red,1)
    Ftest.F(i,:) = fdiff(redMod.RSS_red(i,:), RSS, size(x,2)-1, size(x,2), size(x,1));
end
Ftest.F = max(Ftest.F,0);
Ftest.sign = 1-fcdf(Ftest.F, Ftest.df(1), Ftest.df(2));

end
%% vectorial glm
function [param] = glm_vect(x,y)
covmattr = x'*x;
% eigenvalue decomposition of covariance matrix of full training set
[U, S] = eig(covmattr); % U is matrix of eigenvectors; s is diagonal matrix of eigenvalues
eigvals = diag(S);
Usr = U'*(x'*y); % projection of stim on resp
D = diag(1./eigvals);
param = transpose(U*D*Usr);
end