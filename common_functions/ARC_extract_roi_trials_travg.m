function R = ARC_extract_roi_trials_travg(subjDat,ri,cfg,si)
roiMaskFile = cfg.anatMasks{ri};
roiMaskVol  = spm_read_vols(spm_vol(fullfile(subjDat.subjectDir,roiMaskFile)));
roiMaskVol(isnan(roiMaskVol)) = 0;
R.mask = logical(roiMaskVol);

glm_mask = logical(spm_read_vols(spm_vol(fullfile(subjDat.subjectDir,'ARC3_anatgw.nii'))));
anat3 = R.mask(glm_mask);

% ---------- single-trial betas ----------
glmFile = fullfile(subjDat.subjectDir,'full_zsc_sesswise.mat');
glmData = load(glmFile,'full_zsc');
mc = squeeze(glmData.full_zsc(logical(anat3),:));      % [vox x trial]
[r,c] = find(isnan(mc));
mc(r,:) = [];

R.betatrials = mc;
R.behav = subjDat.behav;                       % keep behavioural ratings
R.group_vec = repmat((1:160)',3,1); % Trial wise structure of which odor was delivered at each trial   
R.si = si;
end
