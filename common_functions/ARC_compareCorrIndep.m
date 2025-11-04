function [zobs,p,opts] = ARC_compareCorrIndep(r1,r2,df1,df2,varargin)
% compareCorrIndep  Compare two independent Pearson r's with Fisher z-test
% Inputs:
%   r1,r2    - correlation coefficients (scalars)
%   df1,df2  - degrees of freedom for each correlation (assumed = n-2)
% Optional name/value:
%   'Tail'   - 'two' (default), 'right', or 'left' for one/two-sided test
%
% Outputs:
%   zobs     - z statistic
%   p        - p-value (according to Tail)
%   opts     - struct with intermediate values (n1,n2,z1,z2,SE)

p = inputParser;
addParameter(p,'Tail','two',@(s) ismember(s,{'two','right','left'}));
parse(p,varargin{:});
tail = p.Results.Tail;

% convert df -> n (assumes df = n - 2)
n1 = df1 + 2;
n2 = df2 + 2;
if n1 <= 3 || n2 <= 3
    error('Sample sizes too small (n must be > 3). Check your df inputs.');
end

z1 = atanh(r1);
z2 = atanh(r2);
SE = sqrt( 1/(n1-3) + 1/(n2-3) );

zobs = (z1 - z2) / SE;

switch tail
    case 'two'
        p = 2*(1 - normcdf(abs(zobs)));
    case 'right'
        p = 1 - normcdf(zobs);   % H1: r1 > r2
    case 'left'
        p = normcdf(zobs);       % H1: r1 < r2
end

opts = struct('n1',n1,'n2',n2,'z1',z1,'z2',z2,'SE',SE);
end
