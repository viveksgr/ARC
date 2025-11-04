function [beta_mat, p_vals, info] = ARC_diff2conds_lessConservative(beta_array, dfmat, varargin)
% beta_array: [S x R x 2] Pearson r per subject, ROI, condition
% dfmat     : [S x R x 2] df (= n-2) for each r
% Options:
%   'Rho'    : assumed corr between z1 and z2 within a subject/ROI (default 0.3)
%   'Method' : 'fixed' (default, meta-analytic Z) or 't' (one-sample t on per-subj z)
%
% Outputs:
%   beta_mat : [S x R] per-subject normalized diff  z_sr = (z2 - z1)/sqrt(VarDiff)
%   p_vals   : [1 x R] p-values per ROI (according to Method)
%   info     : struct with z_fixed, p_fixed, t, p_t, dz, VarDiff, weights

% ---- options ----
p = inputParser;
addParameter(p,'Rho',0.3,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'Method','fixed',@(s) ismember(lower(s),{'fixed','t'}));
parse(p,varargin{:});
rho    = p.Results.Rho;
method = lower(p.Results.Method);

[S,R,C] = size(beta_array);
if C~=2 || ~isequal(size(dfmat),[S R C])
    error('beta_array and dfmat must be SxRx2 with matching sizes.');
end

beta_mat = nan(S,R);            % per-subject normalized difference
dz       = nan(S,R);            % Fisher-z differences
VarDiff  = nan(S,R);            % variance of z2 - z1
W        = nan(S,R);            % inverse-variance weights

for s = 1:S
    for r = 1:R
        r1  = beta_array(s,r,1); r2  = beta_array(s,r,2);
        df1 = dfmat(s,r,1);      df2 = dfmat(s,r,2);
        if ~all(isfinite([r1 r2 df1 df2])), continue; end
        n1 = df1 + 2; n2 = df2 + 2; if n1<=3 || n2<=3, continue; end
        z1 = atanh(r1); z2 = atanh(r2);
        v  = 1/(n1-3) + 1/(n2-3) - 2*rho/sqrt((n1-3)*(n2-3)); % dependence-adjusted
        if v<=0 || ~isfinite(v), continue; end
        dz(s,r)        = z2 - z1;
        VarDiff(s,r)   = v;
        W(s,r)         = 1/v;
        beta_mat(s,r)  = dz(s,r) / sqrt(v);                    % per-subject z
    end
end

p_vals  = nan(1,R);
z_fixed = nan(1,R);
t_stat  = nan(1,R);
p_t     = nan(1,R);

for r = 1:R
    mask = isfinite(dz(:,r)) & isfinite(W(:,r)) & W(:,r)>0;
    if ~any(mask), continue; end

    % Fixed-effect meta-analytic Z on Fisher-z differences (powerful for small S)
    num = sum(W(mask,r) .* dz(mask,r));
    den = sqrt(sum(W(mask,r)));
    z   = num / den;
    z_fixed(r) = z;
    p_fixed    = 2*(1 - normcdf(abs(z)));

    % One-sample t on per-subject normalized z (more conservative with small S)
    zsr = beta_mat(mask,r);
    m   = numel(zsr);
    if m>=2
        mu  = mean(zsr);
        sd  = std(zsr,0);
        t   = mu / (sd / sqrt(m));
        ptt = 2*(1 - tcdf(abs(t), m-1));
        t_stat(r) = t; p_t(r) = ptt;
    end

    p_vals(r) = strcmp(method,'fixed') * p_fixed + strcmp(method,'t') * p_t(r);
end

info = struct('z_fixed',z_fixed,'p_fixed',2*(1-normcdf(abs(z_fixed))), ...
              't',t_stat,'p_t',p_t,'dz',dz,'VarDiff',VarDiff,'weights',W);
end
