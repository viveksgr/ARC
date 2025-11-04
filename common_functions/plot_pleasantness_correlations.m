function out = plot_pleasantness_correlations(behav, featLabels, varargin)
% plot_pleasantness_correlations
%   Bar plot of mean correlations between pleasantness (col 2) and all other
%   features, averaged across subjects (ignoring NaNs), sorted descending.
%   Individual subject correlations are overlaid as points.
%
% Inputs
%   behav       : struct array, one per subject
%                 behav(s).ratings [N x F], behav(s).cid (unused here)
%                 Assumes F >= 3 and pleasantness is column 2 for all subjects.
%   featLabels  : 1 x (F-1) cell array of labels for all features EXCEPT col 2,
%                 ordered by column index (i.e., feature 1,3,4,...,F).
%
% Name/Value (optional)
%   'Title'     : figure title (default 'Correlations with Pleasantness')
%   'MarkerSize': size for subject dots (default 36)
%   'BarColor'  : RGB for bars (default [0.4 0.6 0.85])
%   'DotColors' : S x 3 RGB matrix for subjects (default grayscale)
%
% Output (struct)
%   out.R          : [S x M] subject-by-feature correlations (M = F-1)
%   out.meanR      : [1 x M] across-subject mean correlations (omitnan)
%   out.sortIdx    : permutation applied (descending mean)
%   out.labelsSorted : labels after sorting
%   out.axes       : axes handle
%   out.fig        : figure handle

p = inputParser;
addParameter(p,'Title','Correlations with Pleasantness');
addParameter(p,'MarkerSize',36);
addParameter(p,'BarColor',[0.4 0.6 0.85]);
addParameter(p,'DotColors',[]);
parse(p,varargin{:});
opt = p.Results;

S = numel(behav);
assert(S >= 1, 'behav must have at least one subject.');
F = size(behav(1).ratings, 2);
assert(F >= 3, 'Expected at least 3 features (pleasantness is col 2).');

% Indices of features excluding column 2
otherIdx = setdiff(1:F, 2, 'stable');
M = numel(otherIdx);
assert(numel(featLabels) == M, 'featLabels must be 1x%d for all non-pleasantness features.', M);

% Collect correlations per subject
R = nan(S, M);
for s = 1:S
    X = behav(s).ratings;
    if size(X,2) ~= F
        error('Subject %d has %d columns, expected %d.', s, size(X,2), F);
    end
    y = abs(X(:,2)); % pleasantness
    for j = 1:M
        k = otherIdx(j);
        r = corr(y, X(:,k), 'Rows','pairwise', 'Type','Pearson');
        R(s,j) = r; % may be NaN if insufficient non-NaN overlap
    end
end

% Across-subject mean (ignore NaNs per feature)
meanR = mean(R, 1, 'omitnan');

% Sort descending by mean correlation
[meanR_sorted, sortIdx] = sort(meanR, 'descend');
R_sorted = R(:, sortIdx);
labels_sorted = featLabels(sortIdx);

% Optionally drop features with all-NaN across subjects (nothing to show)
good = ~isnan(meanR_sorted);
meanR_sorted = meanR_sorted(good);
R_sorted = R_sorted(:, good);
labels_sorted = labels_sorted(good);
Mgood = numel(meanR_sorted);

% Colors for subject dots
if isempty(opt.DotColors)
    % soft grayscale per subject
    cmap = linspace(0.2, 0.6, max(S,2))';
    dotCols = [cmap, cmap, cmap];
else
    dotCols = opt.DotColors;
    if size(dotCols,1) < S
        error('DotColors must have at least S=%d rows.', S);
    end
end

% Plot
out.fig = figure('Color','w');
ax = axes('Parent', out.fig); hold(ax,'on');

hb = bar(ax, 1:Mgood, meanR_sorted, 'FaceColor', opt.BarColor, 'EdgeColor','none');
ylim(ax, [-1 1]);
xlim(ax, [0.5, Mgood+0.5]);
grid(ax,'on'); box(ax,'off');

% Overlay subject dots (with tiny jitter)
jitter = linspace(-0.15, 0.15, S);
for s = 1:S
    xs = (1:Mgood) + jitter(s);
    scatter(ax, xs, R_sorted(s, :), opt.MarkerSize, dotCols(s,:), 'filled', ...
        'MarkerEdgeColor','w');
end

% Cosmetics
xticks(ax, 1:Mgood);
xticklabels(ax, labels_sorted);
xtickangle(ax, 40);
ylabel(ax, 'Pearson r');
title(ax, opt.Title);

% Legend
legStr = arrayfun(@(s)sprintf('S%d',s), 1:S, 'uni', 0);
legend(ax, [{'Mean'}, legStr{:}], 'Location','bestoutside');

% Bring bars behind points
uistack(hb,'bottom');

% Output
out.R = R;
out.meanR = meanR;
out.sortIdx = sortIdx(good);
out.labelsSorted = labels_sorted;
out.axes = ax;
end

% 
% nS = 3;
% behav_lab1 = behav(1).percepts;
% behav_lab2 = behav(2).percepts;
% for ii=1:nS
%     behav(ii).rel(end+1)=0;
% end
% 
% idx2 = [1:10 19 11:14 19 15 19 16:18]; % Use this is argsort in sub(2 and 3) to match the labels
% idx1 = [1:18 19 19 19];
% behav_labs = {behav_lab1{:} behav_lab2{end-2:end}};
% 
% % Reordering.
% behav(1).ratings(:,19) = 0;
% behav(1).ratings = behav(1).ratings(:,idx1);
% 
% % Reordering.
% behav(2).ratings(:,19) = 0;
% behav(2).ratings = behav(2).ratings(:,idx2);
% 
% % Reordering.
% behav(3).ratings(:,19) = 0;
% behav(3).ratings = behav(3).ratings(:,idx2);