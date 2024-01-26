function M_new = ARC_binAndTransformQuantiles(M, V, b)
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

    % Bin the values in V
    [~,~,bin] = histcounts(V, edges);

    % Initialize the new matrix
    [n, ~] = size(M);
    M_new = zeros(n, b);

    % Aggregate values in M based on the bins in V
    for i = 1:b
        if any(bin == i)
            M_new(:, i) = mean(M(:, bin == i), 2);
        end
    end
end
