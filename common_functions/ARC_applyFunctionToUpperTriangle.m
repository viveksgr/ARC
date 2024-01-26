function modifiedMatrix = ARC_applyFunctionToUpperTriangle(originalMatrix, regressor_)
    % OriginalMatrix is Square
    % regressor is a column vector

    % Ensure the matrix is square
    [rows, cols] = size(originalMatrix);
    if rows ~= cols
        error('Input must be a square matrix');
    end

    % Extract the upper triangular part of the matrix, excluding the diagonal
    triuIndices = triu(true(size(originalMatrix)), 1);
    upperTriangularElements = originalMatrix(triuIndices);

    % Apply the custom function to these elements
    modifiedElements = regressmeout(upperTriangularElements',regressor_')';

    % Initialize the matrix that will store the modified elements
    modifiedMatrix = originalMatrix;

    % Assign the modified elements back to their original positions
    modifiedMatrix(triuIndices) = modifiedElements;
end
