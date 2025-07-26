function [res, vExp] = ARC_voxelPLS(X, v1_vec, v2_vec)
% ARC_voxelPLS  Supervised PLS to extract trial scores for two behaviours
%
% INPUT
%   X   : N × P  (trials × voxels)  -- raw voxel patterns
%   v1_vec  : N × 1  ordinal / continuous behaviour 1
%   v2_vec  : N × 1  behaviour 2
%
% OUTPUT
%   scores : N × 2   latent scores (one column per behaviour)
%   w      : P × 2   voxel weight vectors
%   vExp   : 1 × 2   variance in Y explained by each component
%
% -------------------------------------------------------------------
% 1) z-score predictors and responses
Xz  = zscore(X);                      % centre & scale columns
Yz  = zscore([v1_vec(:) v2_vec(:)]);          % N × 2
% 2) two-component PLS
nComp = 2;
[c,~,XS,~,beta,PCTVAR] = plsregress(Xz,Yz,nComp);
% 3) trial scores are XS (N × nComp)
scores = XS(:,1:2);
% 4) voxel weights: first nComp rows of beta (skip intercept)
w = beta(2:end,1:2);                  % P × 2

% w = c;
% 5) variance explained in Y
vExp = 100*PCTVAR(2,1:2);             % (%)

res.w_scores = w;
res.t_corr = fastcorr(w(:,1),w(:,2));
end
