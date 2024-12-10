function mdl = ARC_regress_nested_mdl(nmatred, grp, nfolds)
if nargin < 3
    nfolds = 10; % Default to 10-fold cross-validation if not specified
end

nfeatures = size(nmatred,2);
% Define the training set for this fold
trainX = nmatred;
trainY = grp;

% Inner cross-validation for parameter tuning
inner_trainind = crossvalind('Kfold', length(trainY), nfolds); % Same number of folds
bestC = 1; % Default C
bestMSE = inf; % Track the best (lowest) MSE

% Parameter grid for C (broader and more granular)
C_vals = logspace(-2, 1, 7);  % Example: spans from 0.01 to 10 in logarithmic scale
kernel_scales = 1/nfeatures*2 .^ (-6:2:6);
% kernel_scales = 2 .^ (-6:2:6);

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

            % % Train model on inner training set using linear kernel
            % inner_mdl = fitrsvm(inner_trainX, inner_trainY, 'KernelFunction', 'linear', ...
            %                     'BoxConstraint', C, ...
            %                     'Standardize', true);
            %  % Validate model on inner testing set
            % inner_pred = predict(inner_mdl, inner_testX);

            opt = sprintf('-s 3 -t 2 -c %f -g %f -q', C,kernel_scale);
            inner_mdl = svmtrain( inner_trainY,  inner_trainX, opt);
            inner_pred = svmpredict(  inner_testY , inner_testX ,  inner_mdl, ' -q ');

            % Cost function
            mse_inner(inner_idx) = mean((inner_pred - inner_testY).^2);
            % mse_inner(inner_idx) = 1-fastcorr(inner_pred,inner_testY);
        end
        meanMSE = mean(mse_inner);
        if meanMSE < bestMSE
            bestMSE = meanMSE;
            bestC = C;
            bestKernelScale = kernel_scale;
        end
    end
end

% % Train model on the entire training set with best parameters found
% final_mdl = fitrsvm(trainX, trainY, 'KernelFunction', 'linear', ...
%                     'BoxConstraint', bestC, ...
%                     'Standardize', true);
% % Predict on the actual test set
% final_predictions = predict(final_mdl, testX);

opt = sprintf('-s 3 -t 2 -c %f -g %f -q', bestC,bestKernelScale);
mdl = svmtrain( trainY,  trainX, opt);
