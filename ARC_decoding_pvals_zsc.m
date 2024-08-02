function p_values_matrix = ARC_decoding_pvals_zsc(correlations, dfs)
    % Validate input dimensions
    [N, p1, p2] = size(correlations);
    if ~isequal(size(correlations), size(dfs))
        error('The dimensions of correlations and degrees of freedom matrices must match.');
    end
    
    % Initialize the p-values matrix
    p_values_matrix = zeros(p1, p2);
    
    % Reshape the matrices for easier processing
    correlations_reshaped = reshape(correlations, N, p1 * p2);
    dfs_reshaped = reshape(dfs, N, p1 * p2);
    
    % Loop over each element in the p1 x p2 grid
    for idx = 1:(p1 * p2)
        % Extract the vector of correlations and dfs for this grid element
        current_correlations = correlations_reshaped(:, idx);
        current_dfs = dfs_reshaped(:, idx);
        
        % Call the weightedZTest function
        p_values_matrix(idx) = weightedZTest(current_correlations, current_dfs);
    end
    
    % Reshape the p-values back into the p1 x p2 matrix
    p_values_matrix = reshape(p_values_matrix, p1, p2);
end

function p_value = weightedZTest(correlations, dfs)
    % Validate input sizes
    if length(correlations) ~= length(dfs)
        error('The length of correlations and degrees of freedom must be the same.');
    end

    % Step 1: Convert correlation coefficients to Fisher's z-scores
    z_scores = atanh(correlations);  % Fisher's z-transform
    
    % Step 2: Calculate the standard error of each z-score
    ses = 1 ./ sqrt(dfs - 3);  % Standard error of z-scores

    % Step 3: Calculate weights based on the inverse of the variance
    weights = 1 ./ (ses .^ 2);
    
    % Step 4: Calculate the weighted sum of z-scores
    % Including the sign of each z-score in the sum
    weighted_z_scores = z_scores .* weights;
    sum_weighted_z_scores = sum(weighted_z_scores);

    % Step 5: Calculate the variance of the weighted sum
    variance_weighted_z_scores = sum(weights);
    
    % Step 6: Calculate the standard score (z-value) of the weighted sum
    standard_z = sum_weighted_z_scores / sqrt(variance_weighted_z_scores);
    
    % Step 7: Calculate the p-value from the cumulative distribution function
    % for the standard normal distribution (one-tailed test)
    p_value = 1 - normcdf(standard_z);
end

