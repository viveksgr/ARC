function Cavg = ARC_cellAvg3D(Cin)
%ARC_cellAvg3D  Element‐wise average across the 3rd dimension of a cell array
%
%   Cavg = ARC_cellAvg3D(Cin) takes a cell array Cin of size
%   n1×n2×n3, where each Cin{i,j,k} is an m1×m2 numeric matrix of the
%   same size. It returns Cavg of size n1×n2, where
%   Cavg{i,j} = mean( cat(3, Cin{i,j,1},…,Cin{i,j,n3}), 3 ).
%
%   Example:
%     % Create a 2×3×4 cell array, each element is 5×6 random
%     Cin = arrayfun(@(x) {rand(5,6)}, zeros(2,3,4));
%     Cavg = ARC_cellAvg3D(Cin);
%     % Now Cavg is 2×3 cell array; each Cavg{i,j} is 5×6, the mean over the 4.
%
%   See also: cellfun, cat, mean

% Validate input
if ~iscell(Cin)
    error('Input must be a cell array.');
end
[n1,n2,n3] = size(Cin);
Cavg = cell(n1, n2);

for i = 1:n1
    for j = 1:n2
        % Concatenate along 3rd dim
        tmp = cat(3, Cin{i,j,1:n3});
        % Compute element‐wise mean
        Cavg{i,j} = mean(tmp, 3);
    end
end
end
