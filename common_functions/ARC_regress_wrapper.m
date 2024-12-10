function [p_accu, predictions, p_accut] = ARC_regress_wrapper(trainX,trainY,testX,testY, nfolds,nperm)
if nargin<6
    nperm = 1;
end

accumat = zeros(1,nperm);
for pp = 1:nperm

% Normalize training data
mu = mean(trainX);
sigma = std(trainX);
trainX_nn = (trainX - mu) ./ sigma;
% Normalize test data using the same mean and std from training data
testX_nn = (testX - mu) ./ sigma;
mdl = ARC_regress_nested_mdl(trainX_nn,trainY,nfolds);
predictions = svmpredict( testY, testX_nn,  mdl, ' -q ');
accumat(pp)= fastcorr(predictions,testY);
end
p_accu = nanmean(accumat);
[~,p_accut] = ARC_r2t(p_accu,length(testY));
end