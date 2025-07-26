function [M_new,M_new_mat] = ARC_binAndTransform_numctrl_colrow(M, V, b, numdesc,nperm)


if nargin<5
    nperm = 1000;
end

% Bin the values in V
ranger = [min(V) max(V)];
edges = linspace(ranger(1), ranger(2), b+1);
[ncounts,~,bin] = histcounts(V, edges);

if nargin<4
    numdesc = 200;
end

% Initialize the new matrix
[n, ~] = size(M);
M_new_mat = zeros(n, b,nperm);

% Aggregate values in M based on the bins in V

% colmat = zeros(numdesc,nperm);
for pp = 1:nperm
     % perm_idx = datasample(1:length(col_in_bin),numdesc);

    for i = 1:b
        col_in_bin = find(bin == i); % Find columns in the current bin

        if ~isempty(col_in_bin)
            numdesc_op = min(numdesc, length(col_in_bin));
            perm_idx = randperm(length(col_in_bin), numdesc_op);
           
            selected_cols = col_in_bin(perm_idx); % Select random columns
            M_matbin = M(:, selected_cols);
            M_new_mat(:, i,pp) = median(M_matbin, 2);
        end
        colmat()
    end
end

M_new = squeeze(mean(M_new_mat,3));
end