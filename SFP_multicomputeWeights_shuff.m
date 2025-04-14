function [weights_dist] = SFP_multicomputeWeights_shuff(X,y,nperm)

% Remove rows with NaN values in X or y
nanRows = any(isnan([X, y]), 2);
X(nanRows, :) = [];
y(nanRows) = [];

% % Standardize the input variables (X)
X = zscore(X);

weights_dist = zeros(nperm,2);

X = [ones(size(X, 1), 1), X];

% % Standardize the output variable (y)
y = zscore(y);

% % Unnormalized:
% Solve for the weights
for pp = 1:nperm
    y = y(randperm(length(y)));
    weights = X \ y;

    weights_dist(pp,:) = weights(2:3);
end