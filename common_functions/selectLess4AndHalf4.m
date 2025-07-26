function mask = selectLess4AndHalf4(vec)
% SELECTLESS4ANDHALF4   Return a binary mask picking all vec<4 and half of vec==4
%
%   mask = selectLess4AndHalf4(vec) returns a logical vector the same size as
%   vec.  mask(i) is true if vec(i)<4, or if vec(i)==4 and i was randomly
%   selected among half of the “4”s.  If there are an odd number of 4’s, it
%   selects floor(n4/2) of them.
%
% Example:
%   vec  = randi(7,4000,1);
%   mask = selectLess4AndHalf4(vec);
%   sum(vec<4), sum(vec==4), sum(mask)
    % find <4 and ==4
    idx_less4 = vec < 4;
    idx_eq4   = find(vec == 4);

    % how many 4's to pick
    n4 = numel(idx_eq4);
    npick = floor(n4/2);

    % randomly choose half of the 4’s
    sel4 = idx_eq4(randperm(n4, npick));

    % build the mask
    mask = false(size(vec));
    mask(idx_less4) = true;
    mask(sel4)      = true;
end
