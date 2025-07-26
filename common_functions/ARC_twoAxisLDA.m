function [res] = ARC_twoAxisLDA(X, v1, v2)
% ARC_twoAxisLDA  Two separate linear discriminants for binary v1 & v2
%
% INPUT
%   X      N × P   voxel patterns (trials × voxels)
%   v1     N × 1   ordinal / continuous
%   v2     N × 1   ordinal / continuous
%   thr1            (optional) threshold for binarising v1 (default median)
%   thr2            (optional) threshold for binarising v2 (default median)
%
% OUTPUT
%   scores  N × 2   LDA scores per trial  (axis-1, axis-2)
%   voxelW  P × 2   voxel weight maps      (w1 ,  w2)
%   LDAmodels 1 × 2 struct array with fitcdiscr models
%

label1 = v1;
label2 = v2;

% % --- axis 1 -------------------------------------------------------------
% mdl1       = fitcdiscr(X, label1, 'DiscrimType','linear','Delta',0.1);
% coeffs1    = mdl1.Coeffs(1,2);
% w1         = coeffs1.Linear;
% b1         = coeffs1.Const;
% score1     = X*w1 + b1;

deltaList = logspace(-3, 0, 20);
[~, ~, w1] = cvLdaLasso(X, label1, deltaList, 5);

% --- axis 2 -------------------------------------------------------------
% mdl2       = fitcdiscr(X, label2, 'DiscrimType','linear','Delta',0.1);
% coeffs2    = mdl2.Coeffs(1,2);
% w2         = coeffs2.Linear;
% b2         = coeffs2.Const;
% score2     = X*w2 + b2;
[~, ~, w2] = cvLdaLasso(X, label2, deltaList, 5);
% package outputs
% scores      = [score1 score2];    % N × 2
voxelW      = [w1     w2];        % P × 2
% LDAmodels   = [mdl1   mdl2];


res.w_scores = voxelW;
res.t_corr = fastcorr(voxelW (:,1),voxelW(:,2));
end
