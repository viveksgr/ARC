function [p_vals, info] = binop_pvals_by_tail(beta_mat, N_eff, varargin)
% binop_pvals_by_tail
%   Binomial p-values per subject × ROI × tail from % significant voxels.
%
% Inputs
%   beta_mat : [S x R x 2]  percentage OR fraction significant voxels
%                            (:,:,1)=positive, (:,:,2)=negative
%   N_eff    : [S x R]      effective #independent tests per subject×ROI
%
% Name/Value (optional)
%   'AlphaVoxel' : per-voxel alpha used to flag significance (default 0.05)
%   'TwoSidedVoxel' : true if voxelwise test was two-sided (default true)
%   'Alternative' : 'greater' (default), 'less', or 'two-sided'
%   'IsPercent'  : true/false/[](auto). If [], auto-detect (>1 -> percent)
%
% Outputs
%   p_vals  : [S x R x 2] p-values (binomial, per tail)
%   info    : struct with fields:
%             .p0  [1x2] null hit-rate(s) per tail
%             .K   [S x R x 2] tested counts
%             .N   [S x R]     tested N per cell
%             .fraction [S x R x 2] fractions used in test

% ---- options ----
ip = inputParser;
addParameter(ip,'AlphaVoxel',0.05,@(x)isnumeric(x)&&isscalar(x)&&x>0&&x<1);
addParameter(ip,'TwoSidedVoxel',true,@islogical);
addParameter(ip,'Alternative','greater',@(s)ismember(lower(s),{'greater','less','two-sided'}));
addParameter(ip,'IsPercent',[],@(x)islogical(x)||isempty(x));
parse(ip,varargin{:});
alphaV   = ip.Results.AlphaVoxel;
twoSided = ip.Results.TwoSidedVoxel;
alt      = lower(ip.Results.Alternative);
isPct    = ip.Results.IsPercent;

% ---- checks ----
if ndims(beta_mat)~=3 || size(beta_mat,3)~=2
    error('beta_mat must be S x R x 2 (pos,neg).');
end
if ~isequal(size(N_eff), size(beta_mat(:,:,1)))
    error('N_eff must be S x R to match beta_mat(:,:,1).');
end

[S,R,~] = size(beta_mat);

% Determine if values are percentages
if isempty(isPct)
    isPct = any(beta_mat(:) > 1) || max(beta_mat(:)) > 1;  % conservative
end

% Null hit-rate per tail
p0_tail = alphaV / (twoSided + ~twoSided);  % if twoSided=true -> alpha/2, else alpha
p0 = [p0_tail, p0_tail];                    % [pos, neg] same by default

% Prepare outputs
p_vals   = nan(S,R,2);
K        = nan(S,R,2);
fracUsed = nan(S,R,2);

% Round N_eff to nearest valid integer
N = round(N_eff);
N(N < 0) = 0;

% ---- compute p-values per cell ----
for t = 1:2  % 1=positive, 2=negative
    frac = double(beta_mat(:,:,t));
    if isPct, frac = frac / 100; end
    frac(frac < 0) = 0; frac(frac > 1) = 1;

    % Convert to counts with rounding, clamp to [0,N]
    k = round(frac .* N);
    k = max(0, min(k, N));

    % Binomial p-values
    switch alt
        case 'greater'   % Pr[X >= k]
            pv = ones(S,R);
            mask = (N > 0) & isfinite(k);
            pv(mask) = 1 - binocdf(k(mask)-1, N(mask), p0(t));
        case 'less'      % Pr[X <= k]
            pv = ones(S,R);
            mask = (N > 0) & isfinite(k);
            pv(mask) = binocdf(k(mask), N(mask), p0(t));
        case 'two-sided' % 2*min{Pr[X <= k], Pr[X >= k]}
            pv = ones(S,R);
            mask = (N > 0) & isfinite(k);
            if any(mask(:))
                pv_lo = binocdf(k(mask),   N(mask), p0(t));
                pv_hi = 1 - binocdf(k(mask)-1, N(mask), p0(t));
                pv(mask) = min(1, 2*min(pv_lo, pv_hi));
            end
    end

    p_vals(:,:,t)  = pv;
    K(:,:,t)       = k;
    fracUsed(:,:,t)= frac;
end

% Info
info = struct();
info.p0       = p0;
info.K        = K;
info.N        = N;
info.fraction = fracUsed;
info.options  = struct('AlphaVoxel',alphaV,'TwoSidedVoxel',twoSided,...
                       'Alternative',alt,'IsPercent',isPct);
end
