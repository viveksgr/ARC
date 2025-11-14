function p = r2p_onetail(rvals, dfs)
% r2p_onetail
%   Compute one-tailed p-values from Pearson's r and degrees of freedom.
%
% INPUTS
%   rvals : [S x R x C] array of Pearson correlation coefficients
%   dfs   : [S x R x C] array of degrees of freedom (n - 2)
%
% OUTPUT
%   p     : [S x R x C] one-tailed p-values
%
% Formula:
%   t = r * sqrt(df / (1 - r^2))
%   p = 1 - tcdf(t, df)

    if ~isequal(size(rvals), size(dfs))
        error('rvals and dfs must have the same size.');
    end

    % Clamp r to the valid range to avoid Inf t-scores
    rvals = max(min(rvals, 0.999999), -0.999999);

    % Convert r to t
    tvals = rvals .* sqrt(dfs ./ max(1e-12, 1 - rvals.^2));

    % Compute one-tailed p-values
   p = 2 * (1 - tcdf(abs(tvals), dfs));

    % Replace non-finite results with NaN
    p(~isfinite(p)) = NaN;
end
