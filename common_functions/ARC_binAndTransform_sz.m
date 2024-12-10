function M_corr = ARC_binAndTransform_sz(M, V, b, range,nvox)
% Bin the values in V
edges = linspace(range(1), range(2), b+1);
[~,~,bin] = histcounts(V, edges);

% Initialize the new matrix
[n, ~] = size(M);
nperm = 1000;
M_corr_mat = zeros(b, b,nperm);

% Aggregate values in M based on the bins in V
for pp = 1:nperm
    mmat = zeros(nvox,b);
    perm_idx = datasample(1:n, nvox);
    for i = 1:b
       mmat(:, i) = mean(M(perm_idx, bin == i), 2);
    end
      modelmd_binned = zscore(mmat,[],2);
       M_corr_mat (:,:,pp) = corrcoef(modelmd_binned);
end
M_corr = squeeze(mean( M_corr_mat,3));
end