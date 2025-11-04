function [beta_mat, p_vals, debug] = diff2conds_dfNorm_t(beta_array, dfmat, rho)
% diff2conds_dfNorm_t
% Inputs:
%   beta_array : [S x R x 2]  Pearson r per subject, ROI, condition
%   dfmat      : [S x R x 2]  degrees of freedom (df = n-2)
%   rho        : optional assumed corr between the two Fisher-z's
%                (same-subject dependence). Default 0.
%
% Outputs:
%   beta_mat : [S x R]  per-subject normalized difference
%              z_sr = (atanh(r2)-atanh(r1)) / sqrt(VarDiff)
%   p_vals   : [1 x R]  two-sided p-values (one-sample t across subjects)
%   debug    : struct with fields z_diff, VarDiff, n_eff, t, dof, p_t,
%              z_stouffer, p_stouffer (alternate meta-analytic normal test)

if nargin < 3 || isempty(rho), rho = 0; end

[S,R,C] = size(beta_array);
if C ~= 2, error('Expected size(beta_array,3)==2.'); end
if ~isequal(size(dfmat), [S R C]), error('dfmat must match beta_array.'); end

% Pre-allocate
z_diff  = nan(S,R);
VarDiff = nan(S,R);
beta_mat = nan(S,R);

for s = 1:S
    for ridx = 1:R
        r1  = beta_array(s,ridx,1);
        r2  = beta_array(s,ridx,2);
        df1 = dfmat(s,ridx,1);
        df2 = dfmat(s,ridx,2);
        if ~isfinite(r1) || ~isfinite(r2) || ~isfinite(df1) || ~isfinite(df2)
            continue;
        end
        n1 = df1 + 2; n2 = df2 + 2;
        if n1 <= 3 || n2 <= 3, continue; end

        z1 = atanh(r1); z2 = atanh(r2);
        % Variance of difference with optional dependence adjustment:
        % Var(z2 - z1) = 1/(n2-3) + 1/(n1-3) - 2*rho/sqrt((n1-3)*(n2-3))
        v  = 1/(n1-3) + 1/(n2-3) - 2*rho/sqrt((n1-3)*(n2-3));
        if ~isfinite(v) || v <= 0, continue; end

        z_diff(s,ridx)  = z2 - z1;
        VarDiff(s,ridx) = v;
        beta_mat(s,ridx)= (z2 - z1)/sqrt(v);   % per-subject normalized diff
    end
end

% Group test per ROI: one-sample t-test on normalized differences
p_vals = nan(1,R);
t_stat = nan(1,R);
dof    = nan(1,R);
for ridx = 1:R
    zsr = beta_mat(:,ridx);
    v   = isfinite(zsr);
    m   = sum(v);
    if m < 2, continue; end
    mu  = mean(zsr(v));
    sd  = std(zsr(v), 0);                 % sample SD (N-1)
    t   = mu / (sd / sqrt(m));
    p   = 2 * (1 - tcdf(abs(t), m-1));

    p_vals(ridx) = p;
    t_stat(ridx) = t;
    dof(ridx)    = m - 1;
end

% Also provide Stouffer’s Z across subjects as a check (equal weights)
z_st  = nan(1,R); p_st = nan(1,R);
for ridx = 1:R
    zsr = beta_mat(:,ridx);
    v   = isfinite(zsr);
    m   = sum(v);
    if m < 1, continue; end
    z   = sum(zsr(v)) / sqrt(m);
    z_st(ridx) = z;
    p_st(ridx) = 2*(1 - normcdf(abs(z)));
end

debug = struct('z_diff',z_diff, 'VarDiff',VarDiff, ...
               't',t_stat, 'dof',dof, 'p_t',p_vals, ...
               'z_stouffer',z_st, 'p_stouffer',p_st);
end
