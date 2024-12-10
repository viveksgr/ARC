function p_values_3d = ARC_combinePValues3D(t_scores, weights, df)
    % Check if input dimensions match
    if ~isequal(size(t_scores), size(weights))
        error('t_scores and weights must have the same dimensions.');
    end
    
    % Reshape the 3D arrays into 2D
    [N, p1, p2] = size(t_scores);
    t_scores_reshaped = reshape(t_scores, N, p1 * p2);
    weights_reshaped = reshape(weights, N, p1 * p2);
    
    dfs_reshaped = reshape(df,N,p1*p2);
    % Initialize variables
    P = p1 * p2; % Total number of predictors after reshaping
    combined_p_values = zeros(1, P);
    
    % Calculate p-values from individual t-scores
    individual_p_values = 1-tcdf((t_scores_reshaped), dfs_reshaped); % Two-tailed p-values

    % Combine p-values using Fisher's method
    for j = 1:P
        % Sum of -2 * log(p) for each predictor across N subjects
        chi2_stat = -2 * sum(log(individual_p_values(:, j)));
        
        % Degrees of freedom is 2N because -2*log(p) follows a chi-squared distribution with 2 degrees of freedom
        combined_p_values(j) = 1 - chi2cdf(chi2_stat, 2 * N);
    end
    
    % Reshape combined p-values back into the original 3D shape (p1 x p2)
    p_values_3d = reshape(combined_p_values, [p1, p2]);
end
