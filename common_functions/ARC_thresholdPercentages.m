function M = ARC_thresholdPercentages(A, thrs)
%ARC_thresholdPercentages  Compute percent above/below thresholds in cell array
%
%   M = ARC_thresholdPercentages(A, thrs) takes:
%     • A    — an n1×n2 cell array, each A{i,j} a numeric vector
%     • thrs — an n1×n2×2 numeric array of thresholds
%              thrs(:,:,1) are “positive” thresholds
%              thrs(:,:,2) are “negative” thresholds
%
%   It returns:
%     • M    — an n1×n2×2 array, where
%         M(i,j,1) = 100 * (# of elements in A{i,j} > thrs(i,j,1)) / length(A{i,j})
%         M(i,j,2) = 100 * (# of elements in A{i,j} < thrs(i,j,2)) / length(A{i,j})
%
%   Example:
%     A = arrayfun(@(x) {randn(50,1)}, zeros(3,4));      % 3×4 cell of 50×1 vectors
%     th = rand(3,4,2)*2 - 1;                            % random thresholds
%     M = ARC_thresholdPercentages(A, th);
%

[n1, n2] = size(A);
assert(isequal(size(thrs), [n1, n2, 2]), ...
    'Threshold array must be size %dx%dx2.', n1, n2);

% Preallocate result
M = zeros(n1, n2, 2);

% Loop over cells
for i = 1:n1
    for j = 1:n2
        vec = A{i,j}(:);
        N   = numel(vec);
        if N==0
            pctAbove = NaN;
            pctBelow = NaN;
        else
            % positive threshold
            thrPos = thrs(i,j,1);
            pctAbove = 100 * sum(vec > thrPos) / N;
            % negative threshold
            thrNeg = thrs(i,j,2);
            pctBelow = 100 * sum(-vec > thrNeg) / N;
        end
        M(i,j,1) = pctAbove;
        M(i,j,2) = pctBelow;
    end
end
end
