function h = plot_varset_roi_lines(var_set, varargin)
% var_set: [S x R x C x P]  (subjects × ROIs × conditions(=2) × parameters)
% Optional name/value:
%   'RoiLabels'   - 1xR cellstr
%   'CondLabels'  - 1x2 cellstr (default {'Cond 1','Cond 2'})
%   'ParamLabels' - 1xP cellstr
%   'YLabel'      - char
%   'Title'       - char

p = inputParser;
addParameter(p,'RoiLabels',[]);
addParameter(p,'CondLabels',{'Cond 1','Cond 2'});
addParameter(p,'ParamLabels',[]);
addParameter(p,'YLabel','');
addParameter(p,'Title','');
parse(p, varargin{:});
opt = p.Results;

[S,R,C,P] = size(var_set);
if C ~= 2, error('Expected two conditions (size(var_set,3) == 2).'); end

% Mean across subjects (omit NaNs): R x C x P
M = squeeze(mean(var_set, 1, 'omitnan'));

% Colors
gold   = [0.85 0.66 0.00];   % golden yellow
purple = [0.45 0.20 0.60];   % purple

% Layout
if R == 4, nrows = 2; ncols = 2; else, nrows = 1; ncols = R; end

h.fig = figure('Color','w');
h.ax  = gobjects(R,1);
h.lines = cell(R,2);

for r = 1:R
    h.ax(r) = subplot(nrows, ncols, r); hold on; box off;
    x = 1:P;
    y1 = squeeze(M(r,1,:))';
    y2 = squeeze(M(r,2,:))';

    h.lines{r,1} = plot(x, y1, '-o', 'Color', gold,   'LineWidth', 2, ...
                        'MarkerFaceColor', gold, 'MarkerEdgeColor', 'w');
    h.lines{r,2} = plot(x, y2, '-o', 'Color', purple, 'LineWidth', 2, ...
                        'MarkerFaceColor', purple, 'MarkerEdgeColor', 'w');

    xlim([1 P]); xticks(1:P);
    if isempty(opt.ParamLabels), xticklabels(string(1:P));
    else, xticklabels(opt.ParamLabels); end

    ylabel(opt.YLabel);
    if isempty(opt.RoiLabels), title(sprintf('ROI %d', r));
    else, title(opt.RoiLabels{r}); end
    if r == 1, legend(opt.CondLabels, 'Location','best'); end
    ylim([0 0.2])
end

try, sgtitle(opt.Title); catch, title(h.ax(1), opt.Title); end
end
