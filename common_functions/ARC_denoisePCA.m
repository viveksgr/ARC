function [nmatred_reduced, explainedVariance] = ARC_denoisePCA(nmatred, varianceThreshold)
    if nargin < 2
        varianceThreshold = 90; % Default to 90% variance if not specified
    end
    
    % Perform PCA on the data matrix
    [coeff, score, latent, tsquared, explained] = pca(nmatred);
    
    % Determine the number of components to retain 90% variance
    cumulativeVariance = cumsum(explained);
    numComponents = find(cumulativeVariance >= varianceThreshold, 1, 'first');
    
    % Reduce the data matrix to the number of components that capture the desired variance
    nmatred_reduced = score(:, 1:numComponents);
    
    % Return the explained variance for each of the retained components
    explainedVariance = explained(1:numComponents);
end
