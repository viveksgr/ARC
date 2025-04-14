function idx_balanced = ARC_median_balance_zero(x)
%MEDIAN_BALANCE_ZERO Returns indices of a subset of x with median exactly zero
%
% INPUT:
%   x - a vector of values (should be roughly centered around 0)
%
% OUTPUT:
%   idx_balanced - indices of the subset of x with median exactly 0

% Ensure x is a column vector
x = x(:);

% Find indices of positive, negative, and zero values
pos_idx = find(x > 0);
neg_idx = find(x < 0);
zero_idx = find(x == 0);

n_pos = numel(pos_idx);
n_neg = numel(neg_idx);

% Determine which side is in excess and trim randomly
if n_pos > n_neg
    keep_pos = randsample(pos_idx, n_neg);
    keep_neg = neg_idx;
elseif n_neg > n_pos
    keep_neg = randsample(neg_idx, n_pos);
    keep_pos = pos_idx;
else
    keep_pos = pos_idx;
    keep_neg = neg_idx;
end

% Combine with all zero entries
idx_balanced = sort([keep_pos; keep_neg; zero_idx]);

% Optional: validate
if median(x(idx_balanced)) ~= 0
    warning('Final subset does not have exact median 0 due to numerical precision.');
end

end
