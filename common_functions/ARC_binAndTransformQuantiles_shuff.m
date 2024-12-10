function M_new = ARC_binAndTransformQuantiles_shuff(M, V, b,nshuff)

    % Determine the quantile edges for binning
    edges = quantile(V, linspace(0, 1, b+1));

    % Ensure unique edges for histcounts (add small epsilon to duplicate edges)
    uniqueEdges = unique(edges);
    if length(uniqueEdges) < length(edges)
        for i = 1:length(edges)
            if sum(edges(i) == uniqueEdges) > 1
                edges(i) = edges(i) + eps(edges(i));
            end
        end
    end

    % Initialize the new matrix
    [n, ~] = size(M);

    M_new = zeros(nshuff,n, b);
    for pp = 1:nshuff
        % Bin the values in V
        V = V(randperm(length(V)));
        [~,~,bin] = histcounts(V, edges);
        % Aggregate values in M based on the bins in V
        for i = 1:b
            M_new(pp,:, i) = mean(M(:, bin == i), 2);
        end
    end





end
