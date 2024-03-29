function rsa_val = ARC_RSA_func(modelmd2, absbehav)
    % Compute the VxV correlation matrix for modelmd2
    corrMatrix = corr(modelmd2);
    
    % Compute the VxV similarity matrix based on absbehav
    % This involves creating a matrix where each element (i,j) is the
    % negative absolute difference between elements i and j in absbehav
    V = length(absbehav);
    similarityMatrix = -abs(absbehav - absbehav');
    
    % Extract off-diagonal elements from both matrices
    % Create a logical index for off-diagonal elements
    logicalIndex = ~eye(V, V);
    % offDiagCorr = corrMatrix(logicalIndex);
    % offDiagSim = similarityMatrix(logicalIndex);

       offDiagCorr = extractOffDiagonalBlock(corrMatrix);
    offDiagSim = extractOffDiagonalBlock(similarityMatrix);
    
    % Compute the correlation between off-diagonal entries
    rsa_val = corr(offDiagCorr, offDiagSim);
end


function A2 = extractOffDiagonalBlock(A)
    % Ensure the input is a square matrix and its size is even
    [rows, cols] = size(A);
    if rows ~= cols || mod(rows, 2) ~= 0
        error('Input must be a square matrix with even size.');
    end

    % Calculate the size of the blocks
    blockSize = rows / 2;

    % Extract the off-diagonal block (e.g., top-right block)
    A2 = A(1:blockSize, (blockSize + 1):end);
    A2 = A2(:);

    % Alternatively, to extract the bottom-left block, use:
    % A2 = A((blockSize + 1):end, 1:blockSize);
end
