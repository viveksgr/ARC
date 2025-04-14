function Y_resid = sfp_decorrelate_sig(Y,odorID,reg)
% Suppose you have:
% Y: [N x F] sniff trace data (N trials, F time bins/features)
% odorID: [N x 1], each entry in {1..O} labeling which odor for each trial
% intensity: [O x 1], intensity(odorID) gives the nuisance rating for that odor

Y_resid = NaN(size(Y));
for f = 1:size(Y, 2)  % For each time bin / feature
    y_f = Y(:, f);    % Nx1 data for this feature

    % We want to regress y_f ~ 1 + intensity(odor)
    % Build predictor X: intensity(odorID)
    X = reg(odorID);

    % Fit linear regression: y_f = b0 + b1 * X + e
    %   We can do it quickly with ordinary least squares if we want:
    %   Using backslash:
    A = [ones(size(X)), X];      % Nx2 design matrix
    b = A \ y_f;                  % OLS solution
    e_f = y_f - A*b;             % Residual

    Y_resid(:, f) = e_f;  % Store the residual
end
