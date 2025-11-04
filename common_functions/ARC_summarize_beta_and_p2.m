function [beta_mat, p_vals, r_bar] = ARC_summarize_beta_and_p2(beta, dfmat, agg, aggParam)
% beta:  S×R×C cell, each cell holds vector of r's (voxels/searchlights/etc.)
% dfmat: S×R×C numeric, df for each subject/ROI/cond (df = n-2)
% agg:   aggregator: 'median' | 'fisher_mean' | 'fisher_trim' | 'topq' | 'rms'
% aggParam: struct with fields depending on agg:
%    fisher_trim: .trim = 0.2   (20% total trim, i.e., 10% each tail) 
%    topq:        .q    = 0.1   (top 10% by |r|, signed mean)
%    rms:         .keepSign = true (use sign(mean(r)))
if nargin<3 || isempty(agg), agg = 'fisher_trim'; end
if nargin<4, aggParam = struct(); end

[S,R,C] = size(beta);

% 1) subject-level robust summaries
beta_mat = nan(S,R,C);
for i = 1:S
    for j = 1:R
        for k = 1:C
            x = beta{i,j,k};
            if isempty(x), continue; end
            x = x(:); x = x(isfinite(x));
            if isempty(x), continue; end
            beta_mat(i,j,k) = agg_r(x, agg, aggParam); % << robust aggregator
        end
    end
end

% 2) across-subject inverse-variance Fisher-z test per ROI/cond
p_vals = nan(R,C); r_bar = nan(R,C);
for j = 1:R
    for k = 1:C
        r  = squeeze(beta_mat(:,j,k));
        df = squeeze(dfmat(:,j,k));
        valid = isfinite(r) & isfinite(df);
        if ~any(valid), continue; end
        n  = df(valid) + 2; 
        w  = max(n - 3, 0);
        zi = atanh(r(valid));
        wi = w;
        if ~any(wi>0), continue; end
        zbar = sum(wi.*zi)/sum(wi);
        SE   = 1/sqrt(sum(wi));
        zobs = zbar/SE;
        p_vals(j,k) = 2*(1 - normcdf(abs(zobs)));
        r_bar(j,k)  = tanh(zbar);
    end
end
end

% ---------- helper ----------
function m = agg_r(x, method, p)
switch lower(method)
    case 'median'
        % Fisher-z median then back-transform (robust to many near-zeros)
        % m = tanh(median(atanh(x)));
        m = median(x);
    case 'fisher_mean'
        % Standard Fisher-z mean (less sensitive to tails than raw mean)
        m = tanh(mean(atanh(x)));
    case 'fisher_trim'
        if ~isfield(p,'trim'), p.trim = 0.2; end % total trim fraction
        z  = atanh(x);
        m  = tanh(trimmean(z, 100*p.trim));
    case 'rms'
        keepSign = true; if isfield(p,'keepSign'), keepSign = p.keepSign; end
        sgn = 1;
        if keepSign && mean(x)~=0, sgn = sign(mean(x)); end
        m = sgn * sqrt(mean(x.^2));  % emphasises sparse strong signals
    otherwise
        error('Unknown aggregator "%s".', method);
end
end
