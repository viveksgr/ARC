function [beta_mat, p_vals, r_bar] = ARC_summarize_beta_and_p(beta, dfmat)
% beta:  S x R x C cell, each cell holds vector/array of Pearson r's
% dfmat: S x R x C numeric, degrees of freedom for each r estimate (df = n-2)
% Outputs:
%   beta_mat: S x R x C numeric mean r per subject/ROI/cond
%   p_vals:   R x C two-sided p-values (weighted Fisher z across subjects)
%   r_bar:    R x C combined correlation (tanh(weighted mean z))

[S,R,C] = size(beta);

% 1) subject-level averages
beta_mat = nan(S,R,C);
for i = 1:S
    for j = 1:R
        for k = 1:C
            x = beta{i,j,k};
            if ~isempty(x)
                beta_mat(i,j,k) = mean(x(:), 'omitnan');
            end
        end
    end
end

% 2) across-subject inference per ROI/condition
p_vals = nan(R,C);
r_bar  = nan(R,C);

for j = 1:R
    for k = 1:C
        r  = squeeze(beta_mat(:,j,k));
        df = squeeze(dfmat(:,j,k));
        valid = isfinite(r) & isfinite(df);
        if ~any(valid), continue; end

        n = df(valid) + 2;                 % sample sizes
        w = max(n - 3, 0);                 % inverse-variance weights
        keep = (w > 0) & isfinite(r(valid));
        if ~any(keep), continue; end

        zi   = atanh(r(valid(keep)));      % Fisher z
        wi   = w(keep);
        zbar = sum(wi .* zi) / sum(wi);    % weighted mean z
        SE   = 1 / sqrt(sum(wi));          % SE of zbar
        zobs = zbar / SE;                  % ~ N(0,1) under H0: r=0
        % two-sided p-value (use normcdf if available; fallback to erfc)
        if exist('normcdf','file')
            p = 2 * (1 - normcdf(abs(zobs)));
        else
            p = erfc(abs(zobs) / sqrt(2));
        end

        p_vals(j,k) = p;
        r_bar(j,k)  = tanh(zbar);          % combined effect size (optional)
    end
end
end
