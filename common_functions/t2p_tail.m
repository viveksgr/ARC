function p = t2p_tail(tvals, dfs)
% t2p_onetail
%   Compute one-tailed p-values from t-scores and degrees of freedom.
%
% INPUTS
%   tvals : [S x R x C] array of t-statistics
%   dfs   : [S x R x C] array of degrees of freedom (same size as tvals)
%
% OUTPUT
%   p     : [S x R x C] one-tailed p-values
%
% Formula: p = 1 - tcdf(tvals, dfs)

    if ~isequal(size(tvals), size(dfs))
        error('tvals and dfs must have the same size.');
    end

    % Compute one-tailed p-values
   p = 2 * (1 - tcdf(abs(tvals), dfs));

    % Handle NaNs and out-of-range cases gracefully
    p(~isfinite(p)) = NaN;
end
