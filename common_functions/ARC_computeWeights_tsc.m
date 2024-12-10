function [weights, t_scores] = ARC_computeWeights_tsc(x1, x2, y)

[r,~] = find(isnan([x1,x2,y]));
x1(r,:)=[];
x2(r,:)=[];
y(r,:) = [];


% corrcoef(y,x2)

% Combine the regressors into a design matrix, including a column for the intercept
X = [ones(length(x1), 1), x1(:), x2(:)];

% Solve for the weights
weights = X \ y(:);

% Calculate residuals
residuals = y(:) - X * weights;

% Estimate the variance of the residuals
residual_variance = sum(residuals.^2) / (length(y) - length(weights));

% Calculate the standard error for each weight
standard_errors = sqrt(residual_variance * diag(inv(X'*X)));

% Calculate the t-scores
t_scores = weights ./ standard_errors;
end

