function [t_corr, t_scores] = ARC_multicomputeWeights_tsc_voxwise(m1, m2, anatmat)

nvox = size(anatmat,1);
t_scores = zeros(nvox,2);

binz = size(anatmat,2);
if mod(binz,2)==0; binzpart1 = binz/2; binzpart2 = binzpart1+1; else; binzpart1 = (binz+1)/2 ; binzpart2 = binzpart1; end

for ii = 1:nvox
    anatmat_thisvox = anatmat(ii,:);
    anatcorr =  1-abs(anatmat_thisvox-anatmat_thisvox');

    utl_mask1 = logical(blkdiag(zeros(binz-binzpart1),ones(binzpart1)));
    utl_mask_blk = flipud(utl_mask1);


    [~,tmp] = ARC_multicomputeWeights_tsc( [m1(utl_mask_blk) m2(utl_mask_blk)],anatcorr(utl_mask_blk));
    t_scores(ii,:) = tmp(2:end);
end

tmp = corrcoef(t_scores(:,1),t_scores(:,2));
r_corr = tmp(2);
[~,t_corr] = ARC_r2t(r_corr,length(t_scores(:,1)));