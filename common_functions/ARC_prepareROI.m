function [neural,posIdx,negIdx,mutIdx,labels] = ...
         ARC_prepareROI(cfg,subjIdx,roiIdx,behavVec)
% Load ROI-specific voxel time-courses and derive voxel populations

root = cfg.root;
roiName = cfg.ROIs.names{roiIdx};
maskFile = cfg.ROIs.masks{roiIdx};

% ---- load single-trial beta matrix (vox × trials) ----------------------
subjDir  = fullfile(root,sprintf('ARC%02d',subjIdx),'single',roiName);
tmp = load(fullfile(subjDir,'TYPED_FITHRF_GLMDENOISE_RR.mat'),'modelmd','noisepool');
X = double(squeeze(tmp.modelmd));         % vox × trials
if cfg.noisePoolCut
    X = X(~tmp.noisepool,:);
end
if cfg.sigVoxCut                                    % intersect with fmask
    fMask = spm_read_vols(spm_vol(fullfile(subjDir,'..',cfg.fMaskFile)));
    fMask = fMask(:) > 0; X = X(fMask,:);
end
if cfg.doZscore, X = zscore(X,0,2); end

% ---- voxel-population indices from RSA weights -------------------------
load(cfg.RSAfile,'w_score_mat');

ws = w_score_mat{subjIdx,roiIdx};          % vox × 2 (val+, val- weights)
ws = squeeze(mean(ws,2));
thr = 0.075;                               % user-defined cutoff
posIdx =  ws(:,1) > thr & ws(:,2) <  thr;  % val+ specific vox
negIdx =  ws(:,2) > thr & ws(:,1) <  thr;  % val- specific vox
mutIdx =  ws(:,1) > thr & ws(:,2) >  thr;  % mutually tuned vox

% ---- behavioural labels per trial -------------------------------------
labels    = discretize(behavVec,cfg.bins);
neural    = X.';                           % trials × vox (transpose for libsvm)

end
