function M_new = ARC_binAndTransform_numctrl(M, V, b, range)
    % Bin the values in V
    edges = linspace(range(1), range(2), b+1);
    [ncounts,~,bin] = histcounts(V, edges);
    numdesc = ceil(min(ncounts)/2);
    % Initialize the new matrix
    [n, ~] = size(M);
    nperm = 1000;
    M_new_mat = zeros(n, b,nperm);
    



    % Aggregate values in M based on the bins in V
    for i = 1:b
        col_in_bin = find(bin == i); % Find columns in the current bin
        for pp = 1:nperm
        if ~isempty(col_in_bin)
            numdesc_op = min(numdesc, length(col_in_bin));
            perm_idx = randperm(length(col_in_bin), numdesc_op);
            selected_cols = col_in_bin(perm_idx); % Select random columns
            M_matbin = M(:, selected_cols);
            M_new_mat(:, i,pp) = mean(M_matbin, 2);
        end
        end
    end

    M_new = squeeze(mean(M_new_mat,3));
end