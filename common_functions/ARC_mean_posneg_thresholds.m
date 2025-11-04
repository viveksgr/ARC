function [pos_mean, neg_mean] = ARC_mean_posneg_thresholds(v1, m)
% mean_posneg_thresholds
%   [pos_mean, neg_mean] = mean_posneg_thresholds(v1, [m1 m2])
%   pos_mean = mean of positive values > m1
%   neg_mean = mean of negative values with |value| > m2  (i.e., value < -m2)
%   NaN if no values meet a criterion. NaNs in v1 are ignored.

    if numel(m) ~= 2
        error('Second input must be a two-element vector [m1 m2].');
    end
    m1 = m(1);
    m2 = abs(m(2));  % ensure nonnegative threshold for absolute-value test

    x = v1(:);
    x = x(isfinite(x));

    pos_vals = x(x > 0 & x > m1);
    neg_vals = x(x < 0 & x < -m2);

    % pos_vals = x(x > 0);
    % neg_vals = x(x < 0);

    pos_mean = sqrt(mean(pos_vals.^2, 'omitnan'));
    if isempty(pos_vals), pos_mean = 0; end

    neg_mean = sqrt(mean(neg_vals.^2, 'omitnan'));
    if isempty(neg_vals), neg_mean = 0; end
end
