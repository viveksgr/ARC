function M_new = ARC_binAndTransform_shuffcoarse(M)

if nargin<5
    nshuff = 1000;
end

% Initialize the new matrix
[n, b] = size(M);

M_new = zeros(nshuff,n, b);
for pp = 1:nshuff
    for nn=1:n
    ids = randperm(b);
    M_new(pp,nn, :) = M(nn,ids);
    end
end
end