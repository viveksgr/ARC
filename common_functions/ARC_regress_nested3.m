function [p_accu, predictions, p_accut, bestC] = ARC_regress_nested3(nmatred, grp, nfolds, nperm)
if nargin < 3
    nfolds = 10; % Default to 10-fold cross-validation if not specified
end

if nargin < 4
    nperm = 1; % Default to 1 permutation if not specified
end

% Initialize accumat for storing results from each permutation
accumat = zeros(nperm, 1);

for pp = 1:nperm
    % Initialize predictions array
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

        % Calculate class weights for the training set
        uniqueClasses = unique(trainY);
        classWeights = 1 ./ histc(trainY, uniqueClasses);
        weightStr = arrayfun(@(u, w) sprintf('-w%d %f ', u, w), uniqueClasses, classWeights, 'UniformOutput', false);
        weightStr = [weightStr{:}];

        % Inner cross-validation for parameter tuning
        inner_trainind = crossvalind('Kfold', length(trainY), nfolds);
        bestC = 1; % Initialize best C
        bestMSE = inf; % Initialize best MSE

        % Parameter grid for C
        C_vals = logspace(-2, 1, 7);  % Range of C values

        for C = C_vals
            mse_inner = zeros(nfolds, 1);
            for inner_idx = 1:nfolds
                inner_test_log = (inner_trainind == inner_idx);
                inner_train_log = ~inner_test_log;

                % Split inner training and validation sets
                inner_trainX = trainX(inner_train_log, :);
                inner_trainY = trainY(inner_train_log);
                inner_testX = trainX(inner_test_log, :);
                inner_testY = trainY(inner_test_log);

                % Train the model on the inner training set
                opt = sprintf('-s 3 -t 0 -c %f %s-q', C, weightStr);
                inner_mdl = svmtrain(inner_trainY, inner_trainX, opt);
                inner_pred = svmpredict(inner_testY, inner_testX, inner_mdl, '-q');

                % Calculate MSE as the cost function
                % mse_inner(inner_idx) = mean((inner_pred - inner_testY).^2);
                mse_inner(inner_idx) = 1-fastcorr(inner_pred,inner_testY);
            end
            meanMSE = mean(mse_inner);
            if meanMSE < bestMSE
                bestMSE = meanMSE;
                bestC = C;
            end
        end

        % Train the model on the entire training set with the best parameters found
        opt = sprintf('-s 3 -t 0 -c %f %s-q', bestC, weightStr);
        final_mdl = svmtrain(trainY, trainX, opt);
        final_predictions = svmpredict(testY, testX, final_mdl, '-q');
        predictions(testind_log) = final_predictions;
    end

    % Compute overall performance
    accumat(pp) = fastcorr(predictions, grp);
end

p_accu = mean(accumat);
[~,p_accut] = ARC_r2t(p_accu,length(grp));
end
