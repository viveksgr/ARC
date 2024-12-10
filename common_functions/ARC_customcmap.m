function cmap = ARC_customcmap(num_colors)
    % Define the colors in normalized RGB
    colors = [43, 42, 39; 19, 114, 186; 221, 140, 118; 226, 226, 213] / 255;

    % Check if num_colors is provided, otherwise default to 64
    if nargin < 1
        num_colors = 64;
    end

    % Number of transitions between colors
    num_transitions = size(colors, 1) - 1;

    % Colors per transition
    colors_per_transition = floor(num_colors / num_transitions);

    % Initialize the colormap matrix
    cmap = zeros(num_colors, 3);

    % Current index in colormap
    current_idx = 1;

    for i = 1:num_transitions
        % Start and end color for the current transition
        start_color = colors(i, :);
        end_color = colors(i + 1, :);

        % Linear interpolation between the start and end colors
        for j = 1:colors_per_transition
            if current_idx <= num_colors
                t = (j - 1) / (colors_per_transition - 1);
                cmap(current_idx, :) = (1 - t) * start_color + t * end_color;
                current_idx = current_idx + 1;
            end
        end
    end

    % Ensure the last color is exactly the last specified color
    if current_idx <= num_colors
        cmap(current_idx:end, :) = repmat(colors(end, :), num_colors - current_idx + 1, 1);
    end
end