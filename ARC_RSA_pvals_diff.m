function p_values_matrix = ARC_RSA_pvals_diff(t_scores, betas, dfs)
    % Validate input sizes
    if ~isequal(size(t_scores), size(betas))
        error('t_scores and betas must have the same dimensions.');
    end
    
    [N, p1, p2] = size(t_scores);  % N is typically 3
    p_values_matrix = zeros(p1, 2);
    
    % Iterate over each predictor and calculate the p-value
    kk = 0;
    for idx = [1,3]
        kk = kk+1;
        for ii = 1:p1
            current_t_scores = t_scores(:, ii, idx);
            current_betas = betas(:,ii, idx);
            % Assuming dfs is either a scalar or a vector length N
            t_score_ref = t_scores(:, ii, idx+1);
            beta_ref = betas(:,ii,idx+1);
             p_values_matrix(ii,kk) = computeGroupPValue(current_betas, current_t_scores, beta_ref, t_score_ref, dfs);
        end
        
    end
    
end


function p_value = computeGroupPValue(current_betas, current_t_scores, beta_ref, t_score_ref, dfs);
    % Compute standard errors from t-scores and beta coefficients
    current_ses = abs(current_betas ./ current_t_scores);  % Ensure no division by zero or negative values
    ref_ses =  abs(beta_ref ./ t_score_ref);

    % Calculate weighted beta
    weights = 1 ./ current_ses.^2;
    weightedBeta1 = sum(current_betas .* weights) / sum(weights);

    % Ref beta
    weights_ref = 1 ./ ref_ses.^2;
    weightedBeta2 = sum(beta_ref .* weights_ref) / sum(weights_ref);
    
    % % t_diff
    % % Calculate the difference between the weights
    weightedBeta_diff = weightedBeta1 - weightedBeta2;
    
    % % Calculate the standard error of the difference
    se = sqrt(1 / sum(weights));
    se2 = sqrt(1 / sum(weights_ref));
    se_diff = sqrt(se^2 + se2^2);
    
    % % Calculate the t-score for the difference
    t_score_diff =  weightedBeta_diff / se_diff;
    
    p_value = 2 * tcdf(-abs(t_score_diff), dfs);  % Two-tailed p-value
end

             