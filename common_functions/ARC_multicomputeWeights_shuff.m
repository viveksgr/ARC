function [weights_dist] = ARC_multicomputeWeights_shuff(X,modelmd_binned_shuff,utl_mask,nshuff)

if nargin < 3
    nshuff = 1000;
end
% % Standardize the input variables (X)
X = zscore(X);
X = [ones(size(X, 1), 1), X];

weights_dist = zeros(nshuff,size(X,2));
for ss = 1:nshuff
    X2 = X;
    modelmd_binned = squeeze(modelmd_binned_shuff(ss,:,:));
    modelmd_binned = zscore(modelmd_binned,[],2);
    modelmd_corrcoef = corrcoef(modelmd_binned);
    y = modelmd_corrcoef(utl_mask);

    % Remove rows with NaN values in X or y
    nanRows = any(isnan([X, y]), 2);
    X2(nanRows, :) = [];
    y(nanRows) = [];

  
    % % Standardize the output variable (y)
    y = zscore(y);


    weights_dist(ss,:) = X2 \ y;

    % weights = X2\y;
    % 
    % % Calculate residuals
    % residuals = y - X2 * weights;
    % 
    % % Estimate the variance of the residuals
    % residual_variance = sum(residuals.^2) / (length(y) - size(X2, 2));
    % 
    % % Calculate the standard error for each weight
    % standard_errors = sqrt(residual_variance * diag(inv(X2' * X2)));
    % 
    % % Calculate the t-scores
    % weights_dist(ss,:) = weights ./ standard_errors;
end

end
