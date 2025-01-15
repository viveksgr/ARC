function [weights, t_scores] = ARC_multicomputeWeights_onesc(X, y)
    % Remove rows with NaN values in X or y
    nanRows = any(isnan([X, y]), 2);
    X(nanRows, :) = [];
    y(nanRows) = [];

    % Initialize the weights and t_scores arrays
    weights = zeros(size(X, 2) + 1, 1); % +1 for the intercept
    t_scores = zeros(size(X, 2) + 1, 1);

    % Standardize the output variable (y)
    y = zscore(y);

    % Loop through each column of X
    for i = 1:size(X, 2)
        % Standardize this column of X
        Xi = zscore(X(:, i));
        Xi = [ones(size(Xi, 1), 1), Xi];  % Add intercept

        % Solve for the weights for this predictor
        weights_i = Xi \ y;

        % Calculate residuals
        residuals = y - Xi * weights_i;

        % Estimate the variance of the residuals
        residual_variance = sum(residuals.^2) / (length(y) - size(Xi, 2));

        % Calculate the standard error for each weight
        standard_errors = sqrt(residual_variance * diag(inv(Xi' * Xi)));

        % Calculate the t-scores
        t_scores_i = weights_i ./ standard_errors;

        % Store results for the current column
        weights(i+1) = weights_i(2);  % Skip the intercept for individual weights
        t_scores(i+1) = t_scores_i(2);
    end

    % Separate intercept calculations
    intercept = [ones(size(X, 1), 1), zeros(size(X, 1), 1)];  % Intercept only model
    weights_intercept = intercept \ y;
    residuals_intercept = y - intercept * weights_intercept;
    residual_variance_intercept = sum(residuals_intercept.^2) / (length(y) - size(intercept, 2));
    standard_errors_intercept = sqrt(residual_variance_intercept * diag(inv(intercept' * intercept)));
    t_scores_intercept = weights_intercept ./ standard_errors_intercept;

    % Store intercept in the first position
    weights(1) = weights_intercept(1);
    t_scores(1) = t_scores_intercept(1);
end
