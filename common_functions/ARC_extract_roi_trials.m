function R = ARC_extract_roi_trials(subjDat,ri,cfg,si)
roiMaskFile = cfg.anatMasks{ri};
roiMaskVol  = spm_read_vols(spm_vol(fullfile(subjDat.subjectDir,roiMaskFile)));
roiMaskVol(isnan(roiMaskVol)) = 0;
R.mask = logical(roiMaskVol);

% ---------- single-trial betas ----------
glmFile = fullfile(subjDat.subjectDir,cfg.anatNames{ri}, ...
                   'TYPED_FITHRF_GLMDENOISE_RR.mat');
glmData = load(glmFile,'modelmd');
R.betatrials = squeeze(glmData.modelmd);      % [vox x trial]
R.behav      = subjDat.behav;                       % keep behavioural ratings
R.group_vec = subjDat.group_vec; % Trial wise structure of which odor was delivered at each trial   
R.si = si;
end
