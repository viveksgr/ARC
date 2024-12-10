function plotBarWithSignificance(Y_mean, Y_std, p_value)
    % Check if the input vectors are of the same length
    if length(Y_mean) ~= length(Y_std) || length(Y_mean) ~= length(p_value)
        error('All input vectors must have the same length');
    end

    % Create bar plot
    barHandle = bar(Y_mean);
    hold on;

    % Add error bars
    numBars = length(Y_mean);
    for k = 1:numBars
        errorbar(k, Y_mean(k), Y_std(k), 'k', 'linestyle', 'none');
    end

    % Check significance and add corresponding annotations
    % for i = 1:numBars
    %     yOffset = 0.05* max(Y_mean); % Offset to ensure the stars are just above the bars
    %     if p_value(i) < 0.001
    %         text(i, Y_mean(i) + sign( Y_mean(i))*(Y_std(i)+yOffset), '***', 'HorizontalAlignment', 'center','color',[1 0 0]);
    %     elseif p_value(i) < 0.01
    %         text(i, Y_mean(i) +  sign( Y_mean(i))*(Y_std(i)+yOffset), '**', 'HorizontalAlignment', 'center','color',[1 0 0]);
    %     elseif p_value(i) < 0.05
    %         text(i, Y_mean(i) + sign( Y_mean(i))*(Y_std(i)+yOffset), '*', 'HorizontalAlignment', 'center','color',[1 0 0]);
    %     end
    % end

    for i = 1:numBars
        yOffset = 0.05 * max(abs(Y_mean + Y_std)); % Offset to ensure the stars are just above the bars
        if Y_mean(i) >= 0
            textY = Y_mean(i) + Y_std(i) + yOffset;  % Adjusted y-coordinate for positive Y_mean
        else
            textY = Y_mean(i) - Y_std(i) - yOffset;  % Adjusted y-coordinate for negative Y_mean
        end
        if p_value(i) < 0.001
            text(i, textY, '***', 'HorizontalAlignment', 'center','color',[1 0 0]);
        elseif p_value(i) < 0.01
            text(i, textY, '**', 'HorizontalAlignment', 'center','color',[1 0 0]);
        elseif p_value(i) < 0.05
            text(i, textY, '*', 'HorizontalAlignment', 'center','color',[1 0 0]);
        end
    end

    % hold off;
end
