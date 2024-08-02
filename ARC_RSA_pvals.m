function p_values_matrix = ARC_RSA_pvals(t_scores, betas, dfs)
    % Validate input sizes
    if ~isequal(size(t_scores), size(betas))
        error('t_scores and betas must have the same dimensions.');
    end
    
    [N, p1, p2] = size(t_scores);  % N is typically 3
    p_values_matrix = zeros(p1, p2);
    
    % Reshape t_scores and betas into 2D for processing
    t_scores_reshaped = reshape(t_scores, N, p1 * p2);
    betas_reshaped = reshape(betas, N, p1 * p2);
    
    % Iterate over each predictor and calculate the p-value
    for idx = 1:(p1 * p2)
        current_t_scores = t_scores_reshaped(:, idx);
        current_betas = betas_reshaped(:, idx);
        % Assuming dfs is either a scalar or a vector length N
        p_values_matrix(idx) = computeGroupPValue(current_betas, current_t_scores, dfs);
    end
    
    % Reshape the result back into p1 x p2 matrix
    p_values_matrix = reshape(p_values_matrix, p1, p2);
end


function p_value = computeGroupPValue(betas, t_scores, dfs)
    % Compute standard errors from t-scores and beta coefficients
    ses = abs(betas ./ t_scores);  % Ensure no division by zero or negative values

    % Calculate weighted beta
    weights = 1 ./ ses.^2;
    weightedBeta = sum(betas .* weights) / sum(weights);
    
    % Calculate combined standard error
    combinedSE = sqrt(1 / sum(weights));
    
    % Calculate t-score
    tScore = weightedBeta / combinedSE;
    
    % Compute p-value (using the minimum df for simplicity)
    min_df = min(dfs);
    p_value = 2 * tcdf(-abs(tScore), min_df);  % Two-tailed p-value
end
