function M2 = ARC_average_chunks(M1)
% ARC_average_chunks - Averages every 1000 rows of M1 into one row of M2
%
% Input:
%   M1 - Nx2 matrix where N is a multiple of 1000
%
% Output:
%   M2 - (N/1000)x2 matrix, each row is the mean of a 1000-row chunk from M1

% Validate input
[N, cols] = size(M1);
if mod(N, 1000) ~= 0 || cols ~= 2
    error('Input must be Nx2 with N a multiple of 1000.');
end

n_chunks = N / 1000;
M2 = zeros(n_chunks, 2);

for i = 1:n_chunks
    idx = (1:1000) + (i-1)*1000;
    M2(i, :) = mean(M1(idx, :), 1);
end
end
