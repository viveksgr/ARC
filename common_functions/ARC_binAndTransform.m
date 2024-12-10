function M_new = ARC_binAndTransform(M, V, b, range)
    % Bin the values in V
    edges = linspace(range(1), range(2), b+1);
    [~,~,bin] = histcounts(V, edges);

    % Initialize the new matrix
    [n, ~] = size(M);
    M_new = zeros(n, b);

    % Aggregate values in M based on the bins in V
    for i = 1:b
        M_new(:, i) = mean(M(:, bin == i), 2);
    end
end