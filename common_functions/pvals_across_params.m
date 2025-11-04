function [p_mat, out] = pvals_across_params(var_set, varargin)
% pvals_across_params
%   Test, per ROI×condition, whether there's an effect across parameters.
%   Uses one-way repeated-measures ANOVA over parameters (within-subject),
%   with subject as the repeated factor. Optionally uses Friedman as fallback.
%
% INPUT
%   var_set : [S x R x C x P]  (subjects × ROIs × conditions(=2) × parameters)
%
% NAME/VALUE OPTIONS
%   'FisherZ'   : true/false (default true). If true, apply atanh() to stabilize
%                  variance for correlation-like inputs (values will be clamped).
%   'Method'    : 'rm' (default) | 'friedman'
%                 'rm' = fitrm+ranova repeated-measures ANOVA
%                 'friedman' = nonparametric RM test
%   'ParamLabels' : 1xP cellstr (optional, for bookkeeping)
%
% OUTPUT
%   p_mat  : [R x C] p-values (two-sided) for parameter effect
%   out    : struct with per-ROI/cond diagnostics:
%            .method, .F, .DF1, .DF2, .chi2, .KendallW, .tables
%
% NOTE
%   Any subject with NaN across any parameter (for a given ROI/cond) is dropped.

% ---- parse options ----
p = inputParser;
addParameter(p,'FisherZ',true,@islogical);
addParameter(p,'Method','rm',@(s) ismember(lower(s),{'rm','friedman'}));
addParameter(p,'ParamLabels',[]);
parse(p, varargin{:});
opt = p.Results;
method = lower(opt.Method);

[S,R,C,P] = size(var_set);
if C < 1 || P < 2
    error('Expect at least 1 condition and >=2 parameters.');
end

% Prepare output
p_mat = nan(R, C);
out.method = method;
out.F      = nan(R, C);
out.DF1    = nan(R, C);
out.DF2    = nan(R, C);
out.chi2   = nan(R, C);
out.KendallW = nan(R, C);
out.tables = cell(R, C);   % ranova / friedman tables per cell

% Fisher-z (optional; clamp to avoid Inf)
X = var_set;
if opt.FisherZ
    X = max(min(X, 1-1e-12), -1+1e-12);
    X = atanh(X);
end

for r = 1:R
    for c = 1:C
        M = squeeze(X(:, r, c, :));     % S x P
        if ndims(M) == 1, M = M(:)'; end
        % Drop subjects with any NaN across parameters
        keep = all(isfinite(M), 2);
        M = M(keep, :);
        S_eff = size(M,1);
        if S_eff < 2
            p_mat(r,c) = NaN;
            continue;
        end

        switch method
            case 'rm' % repeated-measures ANOVA via fitrm+ranova
                % Build table: subjects as rows, parameters as vars
                varNames = strcat("Param", string(1:P));
                T = array2table(M, 'VariableNames', cellstr(varNames));
                T.Subject = (1:S_eff)';  % ID (not used in formula, but kept)

                within = table( (1:P)', 'VariableNames', {'Parameter'} );
                if ~isempty(opt.ParamLabels)
                    within.Parameter = categorical(within.Parameter, 1:P, opt.ParamLabels);
                else
                    within.Parameter = categorical(within.Parameter);
                end

                formula = sprintf('Param1-Param%d ~ 1', P);
                try
                    rm  = fitrm(T, formula, 'WithinDesign', within);
                    tbl = ranova(rm, 'WithinModel', 'Parameter');
                    out.tables{r,c} = tbl;

                    % Extract row for 'Parameter'
                    row = strcmp(tbl.Term, 'Parameter');
                    p   = tbl.pValue(row);
                    p_mat(r,c) = p;

                    % Collect F & dfs (if available)
                    if ismember('F', tbl.Properties.VariableNames)
                        out.F(r,c)   = tbl.F(row);
                    end
                    if ismember('DF', tbl.Properties.VariableNames)
                        out.DF1(r,c) = tbl.DF(row); % may be DF for effect
                    end
                    if ismember('DFDen', tbl.Properties.VariableNames)
                        out.DF2(r,c) = tbl.DFDen(row);
                    end
                catch
                    % If fitrm/ranova fails (e.g., no license), fall back to friedman
                    [p, tbl, stats] = friedman(M, 1, 'off');
                    p_mat(r,c) = p;
                    out.tables{r,c} = tbl;
                    % chi2 and Kendall's W
                    chisq = tbl{2,5}; % typical location (depends on MATLAB version)
                    out.chi2(r,c) = chisq;
                    out.KendallW(r,c) = chisq / (S_eff * (P - 1));
                end

            case 'friedman' % nonparametric RM test
                [p, tbl, stats] = friedman(M, 1, 'off');
                p_mat(r,c) = p;
                out.tables{r,c} = tbl;
                % chi2 and Kendall's W
                chisq = tbl{2,5}; % usually chi-square in 2nd row, 5th col
                out.chi2(r,c) = chisq;
                out.KendallW(r,c) = chisq / (S_eff * (P - 1));
        end
    end
end
end
