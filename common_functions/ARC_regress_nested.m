function [p_accu, predictions,p_accut] = ARC_regress_nested(nmatred, grp, nfolds,nperm)
    if nargin < 3
        nfolds = 4; % Default to 4-fold cross-validation if not specified
    end

    % Initialize variables
    predictions = zeros(size(grp));
    pred_val = zeros(nfolds, 1);
    
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
        
        % Inner cross-validation for parameter tuning
        inner_trainind = crossvalind('Kfold', length(trainY), nfolds); % Same number of folds
        bestC = 1; % Default C
        bestKernelScale = 1; % Default Kernel Scale
        bestMSE = inf; % Track the best (lowest) MSE

        % Parameter grid (simple example)
        C_vals = [0.1, 1, 10];
        kernel_scales = [0.1, 1, 10];

        for C = C_vals
            for kernel_scale = kernel_scales
                mse_inner = zeros(nfolds, 1);
                for inner_idx = 1:nfolds
                    inner_test_log = (inner_trainind == inner_idx);
                    inner_train_log = ~inner_test_log;

                    % Split inner training and validation sets
                    inner_trainX = trainX(inner_train_log, :);
                    inner_trainY = trainY(inner_train_log);
                    inner_testX = trainX(inner_test_log, :);
                    inner_testY = trainY(inner_test_log);

                    % Train model on inner training set
                    inner_mdl = fitrsvm(inner_trainX, inner_trainY, 'KernelFunction', 'linear', ...
                                        'BoxConstraint', C, 'KernelScale', kernel_scale, ...
                                        'Standardize', true);
                    % Validate model on inner testing set
                    inner_pred = predict(inner_mdl, inner_testX);
                    % mse_inner(inner_idx) = mean((inner_pred - inner_testY).^2);
                      mse_inner(inner_idx) = 1-fastcorr(inner_pred,inner_testY);
                end
                meanMSE = mean(mse_inner);
                if meanMSE < bestMSE
                    bestMSE = meanMSE;
                    bestC = C;
                    bestKernelScale = kernel_scale;
                end
            end
        end
        
        % Train model on the entire training set with best parameters found
        final_mdl = fitrsvm(trainX, trainY, 'KernelFunction', 'rbf', ...
                            'BoxConstraint', bestC, 'KernelScale', bestKernelScale, ...
                            'Standardize', true);
        % Predict on the actual test set
        final_predictions = predict(final_mdl, testX);
        predictions(testind_log) = final_predictions;
        % pred_val(trl_idx) = sum(final_predictions == testY) / length(testY);
    end
    
    % Compute overall performance
    p_accu = corr(predictions, grp, 'Type', 'Pearson'); % Pearson correlation as accuracy metric
    [~,p_accut] = ARC_r2t(p_accu,length(grp));
end
