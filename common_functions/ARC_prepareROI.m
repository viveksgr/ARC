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
if cfg.searchl
    posIdx = [];
    negIdx = [];
    mutIdx  = [];
else
    if cfg.splithalf_
        if cfg.prctile
            % load(cfg.RSAfile,'prctile_mat');
            ws = cfg.prctile_mat{subjIdx,roiIdx};
        else
            load(cfg.RSAfile,'results_cell');
            ws_struct = results_cell{cfg.seed}.w_score_matdist;
            ws = ws_struct {subjIdx,roiIdx};
        end
    else
        load(cfg.RSAfile,'w_score_mat');
        ws = w_score_mat{subjIdx,roiIdx};
        ws = squeeze(mean(ws,2));

    end
    % vox × 2 (val+, val- weights)

    thr1 = cfg.thrs_1;
    thr2 = cfg.thrs_2;

    % user-defined cutoff
    posIdx =  ws(:,1) > thr1 & abs(ws(:,2)) <  thr2;  % val+ specific vox
    negIdx =  ws(:,2) > thr1 & abs(ws(:,1)) <  thr2;  % val- specific vox
    mutIdx =  ws(:,1) > thr1 & ws(:,2) >  thr1;  % mutually tuned vox
end
% ---- behavioural labels per trial -------------------------------------
if cfg.disc
    labels    = discretize(behavVec,cfg.bins);
else
    labels = behavVec;
end
neural    = X.';                           % trials × vox (transpose for libsvm)

end
