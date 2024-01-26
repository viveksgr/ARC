function M_new = ARC_binAndTransform_sz(M, V, b, range, numdesc)
    % Bin the values in V
    edges = linspace(range(1), range(2), b+1);
    [~,~,bin] = histcounts(V, edges);

    % Initialize the new matrix
    [n, ~] = size(M);
    M_new = zeros(n, b);

    % Aggregate values in M based on the bins in V
    for i = 1:b
        col_in_bin = find(bin == i); % Find columns in the current bin
        if ~isempty(col_in_bin)
            numdesc_op = min(numdesc, length(col_in_bin));
            perm_idx = randperm(length(col_in_bin), numdesc_op);
            selected_cols = col_in_bin(perm_idx); % Select random columns
            M_matbin = M(:, selected_cols);
            M_new(:, i) = mean(M_matbin, 2);
        end
    end
end