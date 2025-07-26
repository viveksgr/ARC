function ARC_boxplot_groups(dataCell, labels)
%BOXPLOT_GROUPS  Plot grouped boxplots in pairs with gaps between pairs
%
%   boxplot_groups(dataCell) takes a 1×n cell array where each cell
%   contains a numeric vector. It plots n boxplots such that:
%     • boxes 1 & 2 are adjacent,
%     • boxes 3 & 4 are adjacent,
%     • boxes 5 & 6 are adjacent, etc.,
%     with a small gap between each pair.
%
%   boxplot_groups(dataCell, labels) also sets the x-tick labels.
%
% Example:
%   data = {randn(20,1), randn(20,1)+1, randn(20,1)+2, randn(20,1)+3};
%   labels = {'A1','A2','B1','B2'};
%   boxplot_groups(data, labels);

% Validate input
n = numel(dataCell);
if nargin < 2
    labels = cell(1,n);
end
assert(iscell(dataCell) && numel(labels)==n, ...
    'dataCell must be 1×n cell and labels must be 1×n cell.');

% Compute x positions: 1,2, 4,5, 7,8, ...
positions = (1:n) + floor((0:n-1)/2);

% Flatten data and group indices
values = [];
groups = [];
for i = 1:n
    x = dataCell{i}(:);
    values = [values; x];
    groups = [groups; repmat(i, numel(x), 1)];
end

% Create the boxplot
figure;
boxplot(values, groups, ...
    'Positions', positions, ...
    'Widths', 0.6, ...
    'Colors', lines(n), ...
    'Symbol','k+');
hold on;

% Adjust x-ticks and labels
set(gca, 'XTick', positions, 'XTickLabel', labels);
xlim([min(positions)-1, max(positions)+1]);
ylabel('Value');
title('Grouped Boxplots (in pairs)');

% Add vertical separators between pairs
maxY = ylim; 
for g = 1:floor(n/2)-1
    xsep = positions(2*g) + diff(positions(2*g+[0 1]))/2;
    plot([xsep xsep], maxY, 'k--', 'HandleVisibility','off');
end

hold off;
end
