function [rcrit, tcrit] = ARC_p2r(p, df_t)
% p2r_twotailed  Map a two-tailed p-value to the critical Pearson |r|
% 
%   [rcrit, tcrit] = p2r_twotailed(p, df_t)
% 
% Inputs:
%   p     : two-tailed p-value in (0,1)
%   df_t  : t-test degrees of freedom for correlation (df_t = n - 2)
%
% Outputs:
%   rcrit : critical |r| that achieves the given two-tailed p
%   tcrit : corresponding positive t critical value
%
% Note:
%   Sign is symmetric in two-tailed testing. If you need a signed r,
%   just apply sign you want: r = ±rcrit.

    if any(p<0 | p>1)
        error('p must be in (0,1).');
    end
    if any(df_t<=0)
        error('df_t (n-2) must be > 0.');
    end

    % two-tailed: upper-tail probability is p/2
    tcrit = tinv(1 - p/2, df_t);              % positive critical t
    rcrit = tcrit ./ sqrt(tcrit.^2 + df_t);   % |r| corresponding to tcrit
end
