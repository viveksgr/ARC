function [p_accu, predictions,  p_accut, final_mdl] = ARC_regress_normed(nmatred, grp, nfolds, nperm)

if nargin < 3
    nfolds = 10; % Default to 10-fold cross-validation if not specified
end
if nargin < 4
    nperm = 1; % Default to 1 permutation if not specified
end

accumat = zeros(nperm,1);
for pp = 1:nperm
    % Initialize variables
    predictions = zeros(size(grp));

    % Cross-validation setup
    outer_trainind = crossvalind('Kfold', length(grp), nfolds);

    for trl_idx = 1:nfolds
        % Define the test set for this fold
        testind_log = (outer_trainind == trl_idx);
        testX = nmatred(testind_log, :);
        testY = grp(testind_log);

        % Define the training set for this fold
        trainX = nmatred(~testind_log, :);
        trainY = grp(~testind_log);

        % Normalize training data
        mu = mean(trainX);
        sigma = std(trainX);
        trainX_normalized = (trainX - mu) ./ sigma;

        % Normalize test data using the same mean and std from training data
        testX_normalized = (testX - mu) ./ sigma;
       
        opt = sprintf('-s 3 -t 2 -q');
        final_mdl = svmtrain(trainY, trainX_normalized, opt);
        final_predictions = svmpredict(testY, testX_normalized, final_mdl, '-q');

        predictions(testind_log) = final_predictions;
    end

    % Compute overall performance
    accumat(pp) = fastcorr(predictions, grp);
end

p_accu = mean(accumat);
[~, p_accut] = ARC_r2t(p_accu, length(grp));  % Assuming ARC_r2t is a function to compute p-value from correlation
end
