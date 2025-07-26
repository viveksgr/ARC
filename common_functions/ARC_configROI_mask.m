function anatmask = ARC_configROI_mask(cfg,subj,r)

root = cfg.root;
roiName = cfg.ROIs.names{r};
maskFile = cfg.ROIs.masks{r};
subjDir  = fullfile(root,sprintf('ARC%02d',subj),'single');

mask = spm_read_vols(spm_vol(fullfile(subjDir,cfg.maskFile)));
mask(isnan(mask))=0;
mask = logical(mask);

m1 = spm_read_vols(spm_vol(fullfile(subjDir,maskFile)));
m1(isnan(m1)) = 0;
m1 = logical(m1);
anatmask = and( m1,mask);

