function [pct_sig_pos,  thr_pos, pct_sig_neg, thr_neg] = ARC_bidirectionalFDR(r, df, pvals_thr)
% bidirectionalFDR  One‐sided FDR correction for bidirectional effects
%
%   [pct_sig_pos, pct_ns_pos, pct_sig_neg, pct_ns_neg] = ...
%       bidirectionalFDR(r, df, q)
%
% Inputs:
%   r   - Vector of Pearson r values (can be + or –)
%   df  - Degrees of freedom (scalar)
%   q   - Desired FDR level (e.g. 0.05).  Default is 0.05 if omitted.
%
% Outputs:
%   pct_sig_pos - % of all voxels with r>0 that pass the positive‐tail FDR
%   pct_ns_pos  - % of all voxels with r>0 that fail the positive‐tail FDR
%   pct_sig_neg - % of all voxels with r<0 that pass the negative‐tail FDR
%   pct_ns_neg  - % of all voxels with r<0 that fail the negative‐tail FDR
    
    % if nargin<3 || isempty(q), q = 0.05; end
    q = 0.05;
    %--- 1) Compute one‐sided p‐values for each tail
    t = r .* sqrt(df ./ (1 - r.^2));        % convert r to t
    p_pos = 1 - tcdf(t, df);                % P(r > 0)
    p_neg = 1 - tcdf(-t, df);               % P(r < 0)

    %--- 2) Identify which voxels belong to each tail
    pos_idx = (r > 0);
    neg_idx = (r < 0);

    %--- 3) Run BH‐FDR separately on each tail
    if nargin<3
        thr_pos = computeBH( p_pos(pos_idx), q );
        thr_neg = computeBH( p_neg(neg_idx), q );
    else
        thr_pos = pvals_thr(1);
        thr_neg = pvals_thr(2);
    end
    % thr_pos = 0.05;
    % thr_neg = 0.05;

    %--- 4) Compute percentages (relative to *all* voxels)
    N = numel(r);
    pct_sig_pos = sum( pos_idx & (p_pos <= thr_pos) ) / N * 100;
    pct_ns_pos  = sum( pos_idx & (p_pos >  thr_pos) ) / N * 100;
    pct_sig_neg = sum( neg_idx & (p_neg <= thr_neg) ) / N * 100;
    pct_ns_neg  = sum( neg_idx & (p_neg >  thr_neg) ) / N * 100;
    % 
    % %--- (Optional) Display the results
    % fprintf('Positive‐tail FDR (q=%.2f): p ≤ %.4f → %.2f%% voxels significant, %.2f%% non–significant\n',...
    %         q, thr_pos, pct_sig_pos, pct_ns_pos);
    % fprintf('Negative‐tail FDR (q=%.2f): p ≤ %.4f → %.2f%% voxels significant, %.2f%% non–significant\n',...
    %         q, thr_neg, pct_sig_neg, pct_ns_neg);
end

% Helper: simple Benjamini–Hochberg implementation
function thr = computeBH(pvals, q)
    % pvals : vector of p‐values (only the relevant tail)
    % q     : desired FDR level
    if isempty(pvals)
        thr = NaN;
        return
    end
    m = numel(pvals);
    [p_sorted, ~] = sort(pvals);
    crit = (1:m)' * (q/m);
    below = find(p_sorted <= crit);
    if isempty(below)
        thr = 0;      % no p‐values survive
    else
        k   = below(end);
        thr = p_sorted(k);
    end
end
