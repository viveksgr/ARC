function idx = ARC_balanceClasses(labels_val, nbin, samples_per_class)
    % Validate inputs
    if nargin < 3
        samples_per_class = 100; % Default to 100 samples per class if not specified
    end
    
    % Initialize the index vector
    idx = [];

    % Loop through each class
    for k = 1:nbin
        % Find indices of the current class
        current_class_idx = find(labels_val == k);
        
        % Check if current class has enough examples
        num_examples = length(current_class_idx);
        if num_examples == 0
            warning('Class %d has no examples.', k);
            continue;  % Skip this class if no examples are present
        end
        

        sampled_indices = datasample(current_class_idx, samples_per_class);

        % Append the sampled indices to the main index list
        idx = [idx; sampled_indices];
    end
end
