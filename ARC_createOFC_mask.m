fpath = 'C:\Work\ARC\ARC\ARC02\single';
mask1 = fullfile(fpath,'rwofc.nii');
mask2 = fullfile(fpath,'rwvmpfc.nii');
m3 = spm_read_vols(spm_vol(fullfile(fpath,'ARC3_anatgw.nii')));
X1 = load(fullfile(fpath,'OFC','TYPED_FITHRF_GLMDENOISE_RR.mat'),'modelmd');
X1 = squeeze(X1.modelmd);
X2 = load(fullfile(fpath,'VMPFC','TYPED_FITHRF_GLMDENOISE_RR.mat'),'modelmd');
X2 = squeeze(X2.modelmd);
out = combine_masks_features_aligned(mask1, mask2, X1, X2, 'Merge','mean','GMMask',m3);
savepath = fullfile(fpath,'ar_exmats');
mkdir(savepath)
save(fullfile(savepath,'ar_ex_mats.mat'),'out')