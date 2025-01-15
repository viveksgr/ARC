function ARC_plotMatrixLines(M,anat_names)
    % Validate the input matrix dimensions
    if ndims(M) ~= 3 || size(M,2) ~= 4
        error('Input matrix must be n1x4xn3.');
    end

    % Get the size of the matrix
    [n1, ~, n3] = size(M);
    
    % Create a figure to hold the subplots
    figure;
    
    % Loop over each of the 4 slices in the 2nd dimension
    for j = 1:4
        % Create subplot for each slice
        subplot(2, 2, j);
        hold on;
        
        % Get data for this subplot, which is the j-th slice across all n1 and n3
        dataSlice = squeeze(M(:,j,:));
        
        % Check the squeeze operation in case of single dimension anomaly
        if isvector(dataSlice)
            dataSlice = reshape(dataSlice, n1, n3);
        end
        
        % Plot each line within this subplot
        for i = 1:n1
            plot(linspace(-1,1,n3), dataSlice(i, :));
        end
        
        % Add some axis labels and a title
        xlabel('Pleasantness');
        ylabel('RSA (r)');
        title(anat_names{j});
        
        % Optionally set the same scale for all plots or adjust other aesthetics
        axis tight; % Remove excess space
        grid on;
        
        hold off;
    end
    
    % Improve spacing between subplots
    sgtitle('Timeline RSA');
    set(gcf, 'Position', [100, 100, 800, 600]); % Resize figure for better visibility
end
