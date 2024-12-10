% Assuming T1, T2, T3 are already defined as described
% Example modification to include Subject information
T1.Properties.VariableNames = {'CID', 'Odor', 'pls1'};
T2.Properties.VariableNames = {'CID', 'Odor', 'pls2'};
T3.Properties.VariableNames = {'CID', 'Odor', 'pls3'};

% Add missing variables to each table so they can be concatenated
T1.pls2 = NaN(height(T1), 1);
T1.pls3 = NaN(height(T1), 1);
T2.pls1 = NaN(height(T2), 1);
T2.pls3 = NaN(height(T2), 1);
T3.pls1 = NaN(height(T3), 1);
T3.pls2 = NaN(height(T3), 1);

% Combine all tables into one
combinedTable = [T1; T2; T3];

% Group by CID, aggregate Pls values by mean and select any Odor
resultTable = varfun(@nanmean, combinedTable, 'InputVariables', {'pls1', 'pls2', 'pls3'}, ...
    'GroupingVariables', 'CID');

% Since Odor is the same across same CID, we can just grab the first non-NaN entry
odorTable = rmmissing(combinedTable(:, {'CID', 'Odor'}), 'DataVariables', 'Odor');
[~, idx] = unique(odorTable.CID, 'stable');  % Ensure stable unique CID to preserve order
finalOdor = odorTable(idx, :);

% Join tables to combine mean Pls values with corresponding Odor for each CID
finalTable = join(finalOdor, resultTable, 'Keys', 'CID');

% Display the final combined and processed table
disp(finalTable);
