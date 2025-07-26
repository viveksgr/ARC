function blkIdx = blockRowIndices(M)
%BLOCKROWINDICES  Row-to-block map for a block-diagonal matrix.
%
%   blkIdx = BLOCKROWINDICES(M) returns a column-vector whose length equals
%   size(M,1).  blkIdx(i) tells you which diagonal block row *i* of M
%   belongs to.  The first block is numbered 1, the second 2, and so on.
%
%   The routine assumes M is already in block-diagonal form (zeros outside
%   the square blocks along the main diagonal).  Works on full or sparse
%   matrices.

    % ---- basic checks ---------------------------------------------------
    if ndims(M) ~= 2 || size(M,1) ~= size(M,2)
        error('M must be a square 2-D matrix.');
    end

    n = size(M,1);

    % ---- test which rows start new blocks -------------------------------
    % A row starts a new block if everything to the *left* of the diagonal
    % element is zero.  We test that for every row in parallel:
    %
    %   isStart(i) == true  ↔  nnz( M(i,1:i-1) ) == 0
    %
    % Build a logical vector of length n where isStart(i) is true when row i
    % begins a new block.

    % Extract strictly lower-triangular part and count row-wise non-zeros
    lowerCounts = sum(abs(tril(M,-1)) ~= 0, 2);
    isStart     = (lowerCounts == 0);

    % First row always begins the first block (guaranteed true already, but
    % set explicitly for robustness)
    isStart(1)  = true;

    % ---- convert starts to cumulative block numbers ---------------------
    blkIdx      = cumsum(isStart);      % size n-by-1

    % blkIdx(i) now holds the block number of row i
end
