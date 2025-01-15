function ARC_plotColorGradientLines(M)
    % Validate input
    if isempty(M) || ~ismatrix(M)
        error('Input must be a non-empty 2D matrix.');
    end
    
    % Get the size of the matrix
    [n1, n2] = size(M);
    
    % Define the start and end colors (Earthy Yellow to Deep Purple)
    % These colors can be adjusted as per specific RGB values you consider as earthy yellow and deep purple.
    endColor = [247, 202, 24] / 255; % RGB for a bright yellow, adjust as needed
    startColor = [108, 2, 119] / 255; % RGB for a deep purple, adjust as needed

    % Generate color gradient
    colors = [linspace(startColor(1), endColor(1), n2)', ...
              linspace(startColor(2), endColor(2), n2)', ...
              linspace(startColor(3), endColor(3), n2)'];

    % Create a figure
    figure;
    hold on;
    
    % Plot each line with a color from the gradient
    for i = 1:n2
        plot(linspace(-1,1,n1), M(:,i), 'LineWidth', 1, 'Color', colors(i,:));
    end

    % Enhance plot aesthetics
    % title('Integrated domains');
    xlabel('Pleasantness');
    ylabel('Integrated domains');
    hold off;
end
