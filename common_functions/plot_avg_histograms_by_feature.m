function h = plot_avg_histograms_by_feature(behav, featIdx, varargin)
% plot_behav_histograms_avg
%   Make 2x2 histograms of averaged ratings across subjects, aligning by cid.
%
% Inputs
%   behav   : struct array (S=3). behav(s).ratings [Ns x F], behav(s).cid [Ns x 1]
%   featIdx : [3 x 4] matrix. Entry (s,j) = column index in behav(s).ratings to use
%             for subplot j (j=1..4). Use 0 to ignore that subject for that subplot.
%
% Name/Value (optional)
%   'Bins'       : number of bins or vector of bin edges (default 20)
%   'Titles'     : 1x4 cellstr custom subplot titles (default: auto)
%   'XLabel'     : x-axis label for all histograms (default 'Rating')
%   'UseCIDUnion': true/false (default true). If false, uses intersection only.
%
% Output
%   h : struct with figure and axes handles

% ---- params ----
p = inputParser;
addParameter(p,'Bins',20);
addParameter(p,'Titles',[]);
addParameter(p,'XLabel','Rating',@ischar);
addParameter(p,'UseCIDUnion',true,@islogical);
parse(p,varargin{:});
opt = p.Results;

S = numel(behav);
if ~isequal(size(featIdx), [3 4])
    error('featIdx must be 3x4 (subjects x features-to-plot).');
end
if S ~= 3
    warning('Expected 3 subjects; proceeding with S=%d detected.', S);
end

% ---- build union (or intersection) of CIDs across subjects ----
cidLists = cell(S,1);
for s = 1:S
    cidLists{s} = behav(s).cid(:);
end
if opt.UseCIDUnion
    U = unique(vertcat(cidLists{:}));
else
    U = cidLists{1};
    for s=2:S, U = intersect(U, cidLists{s}); end
end
Nu = numel(U);

% ---- pre-align each subject's ratings to U order (by cid) ----
% We’ll build a function that, given a column index k (per subject),
% returns a [Nu x 1] vector aligned to U, with NaN where cid is missing or k==0.
getAlignedCol = @(s,k) aligned_column(behav(s).ratings, behav(s).cid, U, k);

% ---- plotting ----
h.fig = figure('Color','w'); 
h.ax = gobjects(1,4);
titles = opt.Titles;
if isempty(titles), titles = arrayfun(@(j) sprintf('Feature %d', j), 1:4, 'uni', 0); end

for j = 1:4
    % Collect aligned per-subject vectors for this subplot
    Xs = nan(Nu, S);
    for s = 1:S
        k = featIdx(s,j);
        Xs(:,s) = getAlignedCol(s, k);
    end
    % Average across subjects per condition (ignore NaNs, i.e., subjects not contributing)
    Xavg = mean(Xs, 2, 'omitnan');
    Xavg = Xavg(isfinite(Xavg));   % keep only conditions with at least one subject

    % Plot histogram
    h.ax(j) = subplot(2,2,j); 
    histogram(Xavg, opt.Bins, 'FaceColor', [0.3 0.55 0.9]); 
    grid on; box off;
    title(titles{j});
    xlabel(opt.XLabel); ylabel('Count');

    % Small note about how many conditions contributed
    nContrib = numel(Xavg);
    txt = sprintf('N conditions: %d', nContrib);
    xlims = xlim; ylims = ylim;
    text(xlims(1)+0.02*range(xlims), ylims(2)-0.08*range(ylims), txt, ...
         'FontSize', 9, 'Color', [0.2 0.2 0.2]);
end

end

% ===== helper =====
function v = aligned_column(R, cid, U, k)
% R: [Ns x F], cid: [Ns x 1], U: [Nu x 1], k: column index (or 0)
% Output v: [Nu x 1], aligned to U. NaN where subject lacks that cid or k==0.
    Nu = numel(U);
    v = nan(Nu,1);
    if k<=0 || k>size(R,2) || ~isfinite(k), return; end
    [tf, loc] = ismember(U, cid(:));
    idx = find(tf);
    if ~isempty(idx)
        v(idx) = R(loc(tf), k);
    end
end
