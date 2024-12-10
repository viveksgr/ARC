function ARC_compareMatFiles(file1, file2)
    % Load data from both files
    data1 = load(file1);
    data2 = load(file2);

    % Get field names from both structures
    fields1 = fieldnames(data1);
    fields2 = fieldnames(data2);

    % Find common fields and differences
    commonFields = intersect(fields1, fields2);
    uniqueTo1 = setdiff(fields1, fields2);
    uniqueTo2 = setdiff(fields2, fields1);

    % Display unique fields
    if ~isempty(uniqueTo1)
        fprintf('Variables unique to %s:\n', file1);
        disp(uniqueTo1);
    else
        fprintf('There are no unique variables in %s.\n', file1);
    end

    if ~isempty(uniqueTo2)
        fprintf('Variables unique to %s:\n', file2);
        disp(uniqueTo2);
    else
        fprintf('There are no unique variables in %s.\n', file2);
    end

    % Check differences in common fields
    fprintf('\nChecking differences in common variables...\n');
    for i = 1:length(commonFields)
        var1 = data1.(commonFields{i});
        var2 = data2.(commonFields{i});
        
        if isequal(var1, var2)
            continue;  % Skip if variables are identical
        else
            fprintf('Variable %s differs between files.\n', commonFields{i});
            % Optionally, show the values if they are not too large
            if numel(var1) < 10 && numel(var2) < 10  % Arbitrarily chosen size limit for display
                disp('Values in file 1:');
                disp(var1);
                disp('Values in file 2:');
                disp(var2);
            end
        end
    end
end
