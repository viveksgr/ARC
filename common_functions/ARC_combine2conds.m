function [beta_sr, p_vals, r_group] = ARC_combine2conds(beta, dfmat, agg, aggParam)
% combine2conds
%   beta:   S×R×C cell; each cell holds vector of Pearson r's (voxels/etc.)
%   dfmat:  S×R×C numeric; df for each r estimate (df = n-2)
%   agg:    'median' | 'fisher_mean' | 'fisher_trim' | 'topq' | 'rms'  (default: 'fisher_trim')
%   aggParam: struct with fields depending on agg (see agg_r below)
%
% Outputs:
%   beta_sr : S×R numeric, two conditions combined per subject/ROI
%   p_vals  : 1×R two-sided p-values across subjects (H0: r=0)
%   r_group : 1×R inverse-variance–weighted group correlation (optional)

if nargin<3 || isempty(agg), agg = 'fisher_trim'; end
if nargin<4, aggParam = struct(); end

[S,R,C] = size(beta);
if C < 2
    error('Expected two conditions (size(beta,3) >= 2).');
end

% --- 1) subject-level robust summaries per condition ---
r_cell = nan(S,R,C);                 % summarized r per (s,r,c)
for s = 1:S
    for r = 1:R
        for c = 1:C
            x = beta{s,r,c};
            if isempty(x), continue; end
            x = x(:); x = x(isfinite(x));
            if isempty(x), continue; end
            r_cell(s,r,c) = agg_r(x, agg, aggParam);
        end
    end
end

% --- 2) combine conditions per subject/ROI via Fisher-z with weights (n-3) ---
beta_sr = nan(S,R);                  % combined r per (s,r)
w_subj  = zeros(S,R);                % precision weight per (s,r) to use at group level
for s = 1:S
    for r = 1:R
        rc = squeeze(r_cell(s,r,:));
        dfc = squeeze(dfmat(s,r,:));
        valid = isfinite(rc) & isfinite(dfc);
        if ~any(valid), continue; end
        n  = dfc(valid) + 2;
        w  = max(n - 3, 0);          % inverse-variance weights per condition
        if ~any(w>0), continue; end
        zc = atanh(rc(valid));
        z_bar = sum(w .* zc) / sum(w);
        beta_sr(s,r) = tanh(z_bar);
        w_subj(s,r)  = sum(w);       % carry forward total precision for this (s,r)
    end
end

% --- 3) group inference per ROI across subjects (weighted Fisher-z) ---
p_vals  = nan(1,R);
r_group = nan(1,R);
for r = 1:R
    ri = beta_sr(:,r);
    wi = w_subj(:,r);
    valid = isfinite(ri) & (wi > 0);
    if ~any(valid), continue; end
    Zi   = atanh(ri(valid));
    Wi   = wi(valid);
    zbar = sum(Wi .* Zi) / sum(Wi);
    SE   = 1 / sqrt(sum(Wi));
    zobs = zbar / SE;
    if exist('normcdf','file')
        p = 2 * (1 - normcdf(abs(zobs)));
    else
        p = erfc(abs(zobs)/sqrt(2)); % fallback
    end
    p_vals(r)  = p;
    r_group(r) = tanh(zbar);
end
end

% ---------- helper: robust/sparse-aware aggregator ----------
function m = agg_r(x, method, p)
switch lower(method)
    case 'median'
        m = tanh(median(atanh(x)));
    case 'fisher_mean'
        m = tanh(mean(atanh(x)));
    case 'fisher_trim'
        if ~isfield(p,'trim'), p.trim = 0.2; end               % total trim frac
        z  = atanh(x);
        m  = tanh(trimmean(z, 100*p.trim));
    case 'topq'
        if ~isfield(p,'q'), p.q = 0.1; end                      % top 10% |r|
        q  = max(min(p.q,1), eps);
        [~,ord] = sort(abs(x), 'descend');
        k  = max(1, round(q*numel(x)));
        m  = mean(x(ord(1:k)));                                 % signed mean
    case 'rms'
        keepSign = true; if isfield(p,'keepSign'), keepSign = p.keepSign; end
        sgn = 1; mu = mean(x);
        if keepSign && mu~=0, sgn = sign(mu); end
        m = sgn * sqrt(mean(x.^2));
    otherwise
        error('Unknown aggregator "%s".', method);
end
end
