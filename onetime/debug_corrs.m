zz_corrs = zeros(1000,1);

for zz = 1:1000
    modelmd_binneds = squeeze(modelmd_binned_shuff(zz,:,:));
    modelmd_binneds = zscore(modelmd_binneds,[],2);
    modelmd_corrcoef2 = corrcoef(modelmd_binneds);


    modelmd_binned = zscore(modelmd_binned,[],2);


    modelmd_corrcoef = corrcoef(modelmd_binned);
    zz_corrs(zz) = fastcorr(modelmd_corrcoef2(utl_mask2),modelmd_corrcoef(utl_mask2));
end

z = zeros(1000,1); for zz = 1:1000; z(zz) = fastcorr(behav_ratings_,behav_ratings_(randperm(length(behav_ratings_)))); end

z = zeros(100,1); 
for zz = 1:100
    v1 = abs(behav_ratings_-behav_ratings_');
    b2 = behav_ratings_(randperm(length(behav_ratings_)));
    v2 = abs(b2-b2');
       utl_mask2 = logical(triu(ones(length(v2)),1)); % All possible odors
    z(zz) = fastcorr(v1(utl_mask2),v2(utl_mask2)); 
end
