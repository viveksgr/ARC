function voxelLabel = ARC_classifyVoxelSign(neuralMat,alpha,method)
% ARC_classifyVoxelSign
% ---------------------
% Classifies each voxel’s overall response sign across trials.
%
% INPUTS
%   neuralMat :  (N × P) matrix of activity (rows = trials, cols = voxels)
%   alpha     :  significance level for the test (default = 0.05)
%   method    :  't'   – one-sample t-test against zero  (default)
%                'mean'– simple threshold on the mean (no stat test)
%
% OUTPUT
%   voxelLabel : 1 × P integer vector
%                +1  → voxel is significantly positive
%                -1  → voxel is significantly negative
%                 0  → indeterminate (fails significance or |µ|≈0)
%
% Example
%   lbl = ARC_classifyVoxelSign(Y);          % Y = N×P activity matrix
%
% VivekSagar2016, July-2025
% ---------------------------------------------------------------

if nargin < 2 || isempty(alpha),  alpha  = 0.05;  end
if nargin < 3 || isempty(method), method = 't';    end

[N,P]      = size(neuralMat);
voxelLabel = zeros(1,P);           % pre-allocate

switch lower(method)
    case 't'
        % one-sample t-test per voxel (vectorised)
        mu      = nanmean(neuralMat,1);           % 1×P
        sigma   = nanstd (neuralMat,0,1);         % 1×P (N-1 denom)
        se      = sigma ./ sqrt(sum(~isnan(neuralMat),1));
        tval    = mu ./ se;                       % 1×P
        % two-tailed p-value
        pval    = 2 * tcdf(-abs(tval), N-1);
        voxelLabel(pval < alpha &  mu > 0) = +1;  % significantly positive
        voxelLabel(pval < alpha &  mu < 0) = -1;  % significantly negative

    case 'mean'
        thr = 0;                          % simple sign of mean
        mu  = nanmean(neuralMat,1);
        voxelLabel(mu > thr) = +1;
        voxelLabel(mu < thr) = -1;

    otherwise
        error('Unknown method: %s. Use ''t'' or ''mean''.',method);
end
end
