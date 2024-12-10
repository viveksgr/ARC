function p_values_matrix = ARC_decoding_pvals(correlations, dfs)
    % Validate input sizes and dimensions
    [N, p1, p2] = size(correlations);
    if isscalar(dfs)
        dfs = repmat(dfs - 3, N, p1, p2); % Adjust dfs for all elements
    elseif isequal(size(dfs), size(correlations))
        dfs = dfs - 3; % Adjust each element's dfs
    else
        error('dfs must either be a scalar or have the same dimensions as correlations.');
    end

    % Initialize the p-value matrix
    p_values_matrix = zeros(p1, p2);
    
    % Reshape inputs for vectorized processing
    corr_vecs = reshape(correlations, N, p1 * p2);
    dfs_vecs = reshape(dfs, N, p1 * p2);
    
    % Process each vector of correlations
    for idx = 1:(p1 * p2)
        current_corrs = corr_vecs(:, idx);
        current_dfs = dfs_vecs(:, idx);
        
        % Use the existing function to calculate mean correlation and p-value
        [~, p_value] = averageCorrelationsAndTest(current_corrs, current_dfs);
        
        % Store the p-value in the matrix
        p_values_matrix(idx) = p_value; % Fill the p-value matrix
    end
    
    % Reshape the p-value matrix to the correct dimensions
    p_values_matrix = reshape(p_values_matrix, p1, p2);
end

function [mean_r, p_value] = averageCorrelationsAndTest(corr_vec, dfs)
    % Validate input sizes
    if length(corr_vec) ~= length(dfs)
        error('Correlation vector and degrees of freedom vector must be of the same length.');
    end

    % Step 1: Fisher's z-transformation
    z_values = atanh(corr_vec);  % Fisher's z-transform

    % Step 2: Average the Fisher z-values
    mean_z = mean(z_values);

    % Step 3: Transform the average Fisher z-value back to a correlation
    mean_r = tanh(mean_z);

    % % Step 4: Perform a t-test to test if the mean z-value is significantly different from 0
    % % Calculate the standard error of the z-values (SE = 1/sqrt(df-3))
    se_z = sqrt(1 ./ (dfs - 3));

    % % Standard error of the mean z-value
    % sem_z = sqrt(sum(se_z.^2) / length(dfs)^2);
    sem_z = sqrt(1 / sum(1 ./ (se_z.^2)));


    % t-value for the mean z-value
    t_stat = mean_z / sem_z;

    % Degrees of freedom for the combined test
    total_df = sum(dfs - 3);
    % 
    % P-value from t-test
    p_value = 2 * tcdf(-abs(t_stat), total_df);  % Two-tailed p-value
    % p_value = 1 - tcdf(t_stat, total_df);
    % p_value = ARC_r2t(mean_r,mean(dfs));
end

