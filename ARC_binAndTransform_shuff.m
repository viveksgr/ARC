function M_new = ARC_binAndTransform_shuff(M, V, b, range,nshuff)

if nargin<5
    nshuff = 1000;
end

% Initialize the new matrix
[n, ~] = size(M);

M_new = zeros(nshuff,n, b);
for pp = 1:nshuff
    % Bin the values in V
    V = V(randperm(length(V)));
    edges = linspace(range(1), range(2), b+1);
    [~,~,bin] = histcounts(V, edges);
    % Aggregate values in M based on the bins in V
    for i = 1:b
        M_new(pp,:, i) = mean(M(:, bin == i), 2);
    end
end
end