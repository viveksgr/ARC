function [optDelta, mdl, w] = cvLdaLasso(X, labels, deltaList, K)
% cvLdaLasso   Select L1‐penalty for LDA via K‐fold cross‐validation
%
%   INPUTS:
%     X         – [N×P] data matrix
%     labels    – [N×1] binary class labels
%     deltaList – vector of candidate Delta values (e.g. logspace(-3,0,10))
%     K         – number of CV folds (default: 5)
%
%   OUTPUTS:
%     optDelta  – Delta value with lowest CV error
%     mdl       – LDA model fit on all data with optDelta
%     w         – [P×1] weight vector (Coeffs(1,2).Linear)

if nargin<4 || isempty(K)
    K = 5;
end

% Create partition
cvp = cvpartition(labels, 'KFold', K);

losses = zeros(numel(deltaList),1);
for i = 1:numel(deltaList)
    d = deltaList(i);
    foldLoss = zeros(K,1);
    
    for k = 1:K
        trainIdx = training(cvp, k);
        testIdx  = test(cvp, k);
        
        % fit with L1 penalty = d
        mdl_k = fitcdiscr( ...
            X(trainIdx,:), labels(trainIdx), ...
            'DiscrimType','linear', ...
            'Delta', d, ...
            'FillCoeffs','off' ...
        );
        
        yhat = predict(mdl_k, X(testIdx,:));
        foldLoss(k) = mean(yhat ~= labels(testIdx));
    end
    
    losses(i) = mean(foldLoss);
end

% pick best
[~, bestIdx] = min(losses);
optDelta = deltaList(bestIdx);

% retrain on full data
mdl = fitcdiscr( ...
    X, labels, ...
    'DiscrimType','linear', ...
    'Delta', optDelta);

% extract weight vector
coeffs = mdl.Coeffs(1,2);
w = coeffs.Linear;
end
