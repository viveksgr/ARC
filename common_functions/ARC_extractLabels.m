function [labels_val, labels_sal] = ARC_extractLabels(combined_labels, nbins)
    % Extract the original valence and salience labels from the combined labels
    % Assumes combined_labels were calculated as (labels_val - 1) * nbins + labels_sal
    
    % Calculate labels_val and labels_sal from combined_labels
    labels_val = floor((combined_labels - 1) / nbins) + 1;
    labels_sal = mod(combined_labels - 1, nbins) + 1;
end