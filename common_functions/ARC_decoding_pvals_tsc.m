function p_values_matrix = ARC_decoding_pvals_tsc(tscores, dfs)
    % Validate input sizes and dimensions
    [N, p1, p2] = size(tscores);
    if isscalar(dfs)
        dfs = repmat(dfs - 3, N, p1, p2); % Adjust dfs for all elements
    elseif isequal(size(dfs), size(tscores))
        dfs = dfs - 3; % Adjust each element's dfs
    else
        error('dfs must either be a scalar or have the same dimensions as correlations.');
    end

    % Initialize the p-value matrix
    p_values_matrix = zeros(p1, p2);
    
    % Reshape inputs for vectorized processing
    corr_vecs = reshape(tscores, N, p1 * p2);
    dfs_vecs = reshape(dfs, N, p1 * p2);
    
    % Process each vector of correlations
    for idx = 1:(p1 * p2)
        current_corrs = corr_vecs(:, idx);
        current_dfs = dfs_vecs(:, idx);
        
        % Use the existing function to calculate mean correlation and p-value
        [p_value] = ARC_combineDirectionalPValues(current_corrs, current_dfs);
        
        % Store the p-value in the matrix
        p_values_matrix(idx) = p_value; % Fill the p-value matrix
    end
    
    % Reshape the p-value matrix to the correct dimensions
    p_values_matrix = reshape(p_values_matrix, p1, p2);
end





function [combined_p_value] = ARC_combineDirectionalPValues(t_scores, dfs)
    % Calculate z-scores from t-scores and degrees of freedom
    z_scores = t_scores ./ sqrt(dfs);
    
    % Determine the majority direction based on the sign of the z-scores
    majority_direction = sign(sum(sign(z_scores)));
    
    % Adjust z-scores by their agreement with the majority direction
    adjusted_z_scores = z_scores .* (sign(z_scores) == majority_direction);
    
    % Sum the adjusted z-scores
    combined_z_score = sum(adjusted_z_scores);
    
    % Convert the combined z-score to a p-value
    combined_p_value = 1 - normcdf(abs(combined_z_score), 0, 1); % One-tailed p-value
    
    % Check if all tests are in the same direction
    consistent_direction = all(sign(z_scores) == majority_direction);

end
