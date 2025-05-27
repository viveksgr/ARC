function M3 = ARC_inversePercentiles(M1, M2)
%ARC_INVERSEPERCENTILES Compute the percentile rank of means of M1 within M2
%
%   M3 = ARC_inversePercentiles(M1, M2) returns an N×2 matrix M3, where
%   M1 and M2 are both N×1000×2 arrays. For each slice n = 1…N and each
%   variable v ∈ {1,2}, M3(n,v) is the empirical percentile (on [0,1]) of
%   mean(M1(n,:,v)) within the distribution M2(n,:,v).
%
%   Inputs:
%     M1 — N×1000×2 numeric array
%     M2 — N×1000×2 numeric array, same size as M1
%
%   Output:
%     M3 — N×2 array of percentiles (0–1)
%
%   Example:
%     N = 5;
%     M1 = randn(N,1000,2);
%     M2 = randn(N,1000,2);
%     M3 = ARC_inversePercentiles(M1, M2);
%

% Validate inputs
assert(ndims(M1)==3 && ndims(M2)==3, 'M1 and M2 must be 3-D arrays.');
assert(all(size(M1)==size(M2)), 'M1 and M2 must have the same size.');
[N,T,V] = size(M1);
assert(T==1000 && V==2, 'Expected size N×1000×2.');

% Preallocate output
M3 = zeros(N, V);

% Loop over voxels (N) and variables (2)
for n = 1:N
    for v = 1:V
        % Compute the mean of M1 over the 2nd dimension
        m1_mean = mean(M1(n,:,v), 2);  % scalar
        
        % Extract the corresponding bootstrap distribution from M2
        m2_dist = squeeze(M2(n,:,v));  % 1×1000
        
        % Empirical percentile: fraction of M2 entries less than m1_mean
        pct = sum(m2_dist < m1_mean) / numel(m2_dist);
        
        % Assign to output
        M3(n,v) = pct;
    end
end
end
