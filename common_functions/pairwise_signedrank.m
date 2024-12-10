function pairwise_signedrank(data)
    % Input:
    %   data: 4x3 matrix (regions x subjects)
    % Output:
    %   Displays pairwise comparisons with p-values and Bonferroni-adjusted p-values.
    
    % Number of regions
    nRegions = size(data, 1);
    
    % Initialize results
    comparisons = {};
    p_values = [];
    
    % Perform pairwise comparisons
    for i = 1:nRegions
        for j = i+1:nRegions
            % Extract data for the two regions
            region1 = data(i, :);
            region2 = data(j, :);
            
            % Perform Wilcoxon Signed-Rank Test
            p = signrank(region1, region2);
            
            % Store results
            comparisons{end+1} = sprintf('Region %d vs Region %d', i, j);
            p_values(end+1) = p;
        end
    end
    
    % Bonferroni correction
    alpha = 0.05; % Significance level
    nComparisons = length(p_values);
    adjusted_p_values = p_values * nComparisons;
    adjusted_p_values(adjusted_p_values > 1) = 1; % P-values cannot exceed 1
    
    % Display results
    fprintf('\nPairwise Wilcoxon Signed-Rank Test Results:\n');
    fprintf('Comparison\t\t\tp-value\tAdjusted p-value (Bonferroni)\n');
    fprintf('------------------------------------------------------------\n');
    for k = 1:nComparisons
        fprintf('%s\t\t%.4f\t%.4f\n', comparisons{k}, p_values(k), adjusted_p_values(k));
    end
end
