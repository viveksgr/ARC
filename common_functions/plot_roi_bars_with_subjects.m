function h = plot_roi_bars_with_subjects(M, p_vals, varargin)
% M: [num_subj x num_ROI]
% p_vals: [1 x num_ROI] or [num_ROI x 1] (NaN/[] to skip stars)
% Optional name/value:
%   'RoiLabels'   - cellstr, length num_ROI
%   'YLabel'      - char
%   'Title'       - char
%   'BarColor'    - 1x3 rgb (default [0.3 0.55 0.9])
%   'LineColor'   - 1x3 rgb (default [0.2 0.2 0.2])
%   'ErrType'     - 'std' (default) or 'sem'
%   'StarColor'   - 1x3 rgb (default [0 0 0])
%   'StarFontSize'- scalar (default 12)
%   'WeakLevel'   - scalar for '(*)' (default 0.085, set [] to disable)
%   'SigLevels'   - [0.05 0.01 0.001] thresholds
%   'ErrMult'     - multiplier for error bar when placing stars (default 1.5)
%   'StarYOffset' - extra absolute offset added beyond ErrMult*err (default 0)

% ---- parse inputs ----
p = inputParser;
addParameter(p,'RoiLabels',[]);
addParameter(p,'YLabel','');
addParameter(p,'Title','');
addParameter(p,'BarColor',[0.3 0.55 0.9]);
addParameter(p,'LineColor',[0.2 0.2 0.2]);
addParameter(p,'ErrType','std',@(s) ismember(lower(s),{'std','sem'}));
addParameter(p,'StarColor',[0 0 0]);
addParameter(p,'StarFontSize',12);
addParameter(p,'WeakLevel',0.085);
addParameter(p,'SigLevels',[0.05 0.01 0.001]);
addParameter(p,'ErrMult',1.5);
addParameter(p,'StarYOffset',0);
parse(p,varargin{:});
opt = p.Results;

[num_subj, num_roi] = size(M);
if nargin < 2 || isempty(p_vals)
    p_vals = nan(1,num_roi);
else
    p_vals = p_vals(:).'; % make row
    if numel(p_vals) ~= num_roi
        error('p_vals must have length num_ROI.');
    end
end

% stats across subjects
mu = mean(M,1,'omitnan');
sd = std(M,0,1,'omitnan');
switch lower(opt.ErrType)
    case 'std'
        err = sd;
    case 'sem'
        n_eff = sum(isfinite(M),1);
        err = sd ./ max(sqrt(n_eff),1);
end

% ---- plotting ----
figure('Position',[1 1 320 240],'Color','w'); ax = axes; hold(ax,'on');

hb = bar(1:num_roi, mu, 'FaceColor', opt.BarColor, 'EdgeColor','none');
set(hb,'FaceAlpha',0.85);

he = errorbar(1:num_roi, mu, err, 'k', 'LineStyle','none', 'LineWidth',1.2);

% subject lines
for s = 1:num_subj
    plot(1:num_roi, M(s,:), '-o', ...
        'Color', opt.LineColor, 'LineWidth',0.7, ...
        'MarkerSize',4, 'MarkerFaceColor', opt.LineColor, ...
        'MarkerEdgeColor', 'w');
end

% cosmetics
xlim([0.5, num_roi+0.5]); box off;
xticks(1:num_roi);
if isempty(opt.RoiLabels)
    xticklabels(compose('ROI%d',1:num_roi));
else
    xticklabels(opt.RoiLabels);
end
ylabel(opt.YLabel);
title(opt.Title);

% ---- significance stars ----
% decide star string per ROI
sigStr = repmat({''},1,num_roi);
for r = 1:num_roi
    pval = p_vals(r);
    if ~isfinite(pval), continue; end
    if pval < opt.SigLevels(3)
        sigStr{r} = '***';
    elseif pval < opt.SigLevels(2)
        sigStr{r} = '**';
    elseif pval < opt.SigLevels(1)
        sigStr{r} = '*';
    elseif ~isempty(opt.WeakLevel) && pval < opt.WeakLevel
        sigStr{r} = '(*)';
    end
end

% star positions and axis padding
y_text = nan(1,num_roi);
for r = 1:num_roi
    if isempty(sigStr{r}), continue; end
    if mu(r) >= 0
        y_text(r) = mu(r) + opt.ErrMult*err(r) + opt.StarYOffset;
    else
        y_text(r) = mu(r) - opt.ErrMult*err(r) - opt.StarYOffset;
    end
    text(r, y_text(r), sigStr{r}, ...
        'HorizontalAlignment','center', ...
        'VerticalAlignment', ternary(mu(r)>=0,'bottom','top'), ...
        'Color', opt.StarColor, 'FontSize', opt.StarFontSize, ...
        'FontWeight','bold');
end

% ensure stars are visible (expand ylim if needed)
y_all = [mu+err, mu-err, y_text(~isnan(y_text))];
ymin = min(y_all); ymax = max(y_all);
if ~isempty(y_all)
    pad = 0.1 * max(1, ymax - ymin);
    ylim([ymin - pad, ymax + pad]);
end

uistack(he,'top'); % bring error bars to front

% output handles
h = struct('ax',ax,'bar',hb,'err',he);
end

% small helper
function out = ternary(cond,a,b)
if cond, out=a; else, out=b; end
end
