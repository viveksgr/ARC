function ARC_histSplit(cell3vec)
%ARC_histSplit  Plot top-row histograms of 80 shared random samples and
%               bottom-row histograms of the complementary 80 samples.
%
% INPUT
%   cell3vec — 1×3 cell array, each cell a 160×1 numeric vector
%
% The function:
%   • draws the same 80 random indices (without replacement) from 1:160
%   • plots histograms of those samples for the three vectors (top row)
%   • plots histograms of the remaining 80 entries (bottom row)
%
% Example
%   x = randn(160,1); y = randn(160,1)+2; z = randn(160,1)*0.5;
%   ARC_histSplit({x,y,z});

% for ss = 1:3; cell3vec{ss} = behav(ss).ratings(:,2); end

% ── input checks ────────────────────────────────────────────────────────
if ~iscell(cell3vec) || numel(cell3vec)~=3
    error('Input must be a 1×3 cell array of 160×1 vectors.');
end
for k = 1:3
    assert(isvector(cell3vec{k}) && numel(cell3vec{k})==160,...
        'Each cell must contain a 160-element vector.');
    cell3vec{k} = cell3vec{k}(:);               % force column
end

% ── choose indices: 80 in, 80 out ───────────────────────────────────────
idxIn  = sort(randsample(160, 80, false));      % same indices for all
idxOut = setdiff(1:160, idxIn);

% ── create figure & plot ────────────────────────────────────────────────
figure('Color','w'); clf;
titles = {'Vector 1','Vector 2','Vector 3'};

for k = 1:3
    data    = cell3vec{k};
    % top-row : selected 80
    subplot(2,3,k);
    histogram(data(idxIn),7,'FaceColor',[0.2 0.5 0.9]);
    title([titles{k} '  (in 80)']);
    grid on;
    % bottom-row : remaining 80
    subplot(2,3,3+k);
    histogram(data(idxOut),7,'FaceColor',[0.8 0.3 0.3]);
    title([titles{k} '  (out 80)']);
    grid on;
end

sgtitle('Top: random shared 80 samples  |  Bottom: complementary 80');
end
