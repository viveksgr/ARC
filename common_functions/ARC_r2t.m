function [p, t] = ARC_r2t(r, df)
    % Compute the t-statistic from r
    t = (r .* sqrt(df)) ./ sqrt(1 - r.^2);
    
    % Compute the p-value based on the t-distribution
    p = 2 * (1 - tcdf(abs(t), df)); 
end