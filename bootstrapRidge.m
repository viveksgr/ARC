function [beta_m, beta_err, p_values, beta] = bootstrapRidge(X, Y, num_iterations, lambda)
    % X: 160xp predictor matrix
    % Y: 160x1 predicted variable
    % num_iterations: number of bootstrap iterations (e.g., 1000)
    % lambda: ridge regularization parameter

    [num_samples, num_predictors] = size(X);
    
    % Matrix to store beta weights from each bootstrap iteration
    bootstrap_betas = zeros(num_predictors, num_iterations);
    
    for i = 1:num_iterations
        % Sample with replacement
        sample_indices = randsample(num_samples, num_samples, true);
        X_sample = X(sample_indices, :);
        Y_sample = Y(sample_indices);
        
        % Ridge regression
        beta = ridge(Y_sample, X_sample, lambda, 0);  % 0 means no intercept in this function
        bootstrap_betas(:, i) = beta(2:end);  % Exclude the intercept
    end
    
    % Compute mean beta weights across iterations
    beta_m = mean(bootstrap_betas, 2);
    beta_err = std(bootstrap_betas,[],2);
    
    % Compute p-values (using a basic method)
    p_values = zeros(num_predictors, 1);
    for j = 1:num_predictors
        % Count how often beta weight is different from zero
        p_values(j) = computePercentilePValue(bootstrap_betas(j, :)');
    end
end

function p_val = computePercentilePValue(V)

% For a one-tailed test (greater than zero)
p_greater = sum(V > 0) / length(V);

% For a one-tailed test (less than zero)
p_less = sum(V < 0) / length(V);

% For a two-tailed test
p_val = 2 * min(p_greater, p_less);

end
