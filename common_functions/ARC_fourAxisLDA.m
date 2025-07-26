function [res] = ARC_fourAxisLDA(X, v1, v2,v11,v22)
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

% --- axis 1 -------------------------------------------------------------
mdl1       = fitcdiscr(X, label1, 'DiscrimType','linear');
coeffs1    = mdl1.Coeffs(1,2);
w1         = coeffs1.Linear;
b1         = coeffs1.Const;
score1     = X*w1 + b1;

% --- axis 2 -------------------------------------------------------------
mdl2       = fitcdiscr(X, label2, 'DiscrimType','linear');
coeffs2    = mdl2.Coeffs(1,2);
w2         = coeffs2.Linear;
b2         = coeffs2.Const;
score2     = X*w2 + b2;

% axis 3
X_p = X(label1,:);
label11 = v11(label1)';
mdl1       = fitcdiscr(X_p, label11 , 'DiscrimType','linear');
coeffs1    = mdl1.Coeffs(1,2);
w11      = coeffs1.Linear;
% b1         = coeffs1.Const;
% score1     = X*w1 + b1;

% axis 4
X_n = X(label2,:);
label22 = v22(label2)';
mdl1       = fitcdiscr(X_n, label22 , 'DiscrimType','linear');
coeffs1    = mdl1.Coeffs(1,2);
w22      = coeffs1.Linear;

% package outputs
% scores      = [score1 score2];    % N × 2
ab1 = sqrt((zscore(w1).^2+ zscore(w11).^2 )/2);
ab2 = sqrt((zscore(w2).^2+ zscore(w22).^2 )/2);
voxelW      = [ab1     ab2];        % P × 2
% LDAmodels   = [mdl1   mdl2];


res.w_scores = voxelW;
res.t_corr = fastcorr(voxelW (:,1),voxelW(:,2));
end
