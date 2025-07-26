function [Zg,rcrit_pos,rcrit_neg] = ...
           combineRmaps_meta_bi_weighted(rmaps, Ns, q)
% rmaps    : X×Y×Z×S array of each subject’s r-map
% Ns       : S×1 vector of sample sizes (so df_i = Ns(i)-2)
% q        : FDR level (e.g. 0.05)

    if nargin<3, q = 0.05; end
    [X,Y,Z,S] = size(rmaps);
    assert(numel(Ns)==S, 'Need one N per subject');

    % 1) Fisher-Z
    Zmap = atanh(rmaps);

    % 2) Compute weights
    w = sqrt(Ns(:) - 3);        % S×1

    % 3) Voxel-wise weighted sum and normalization
    %    numerator: sum(w_i * Z_i)
    %    denom    : sqrt(sum(w_i^2))
    Wsum = sqrt(sum(w.^2));     % scalar
    Zg = zeros(X,Y,Z);
    for i = 1:S
        Zg = Zg + w(i) * Zmap(:,:,:,i);
    end
    Zg = Zg / Wsum;             % group Z-map

    % write_reshaped_nifty( Zg , pwd, false, fullfile(pwd,'wARC01_val1.nii'),'igmnea.nii');
    Zg(isnan(Zg))=0;
    [~,  thr_pos, ~, thr_neg] = ARC_bidirectionalFDR(Zg, mean(Ns), q);

    rcrit_pos =  r2t(thr_pos,mean(Ns));
    rcrit_neg =  r2t(thr_neg,mean(Ns));
end

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

