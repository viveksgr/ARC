function P = ARC_percentileAbsWeights(W)
% ARC_percentileAbsWeights  Convert weights to column-wise percentiles
%
%   P = ARC_percentileAbsWeights(W)
%
% INPUT
%   W : N × 2   real-valued weight matrix (can contain NaNs)
%
% OUTPUT
%   P : N × 2   percentile ranks (0…100) of |W| within each column
%
% Notes
%   • Uses column-wise tied ranks → uniform percentile definition
%   • NaNs in W propagate to P
%   • Percentile p(i) = 100·(rank-0.5)/n  (Hazen definition)
%
% Author: ChatGPT-o3, Jul 2025
% -------------------------------------------------------------------------
absW = (W);
P    = nan(size(W));

for col = 1:size(absW,2)
    x         = absW(:,col);
    validMask =  ~isnan(x);
    n         =  sum(validMask);
    % tiedrank assigns average rank to ties (1 … n)
    r         = tiedrank(x(validMask));        % column vector
    P(validMask,col) = 100 * (r - 0.5) / n;    % Hazen percentile
end
end
