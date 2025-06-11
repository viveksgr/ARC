function cfg = ARC_makeDecodingCfg(rootDir)
%ARC_makeDecodingCfg  Create a struct with all analysis settings
%
% Usage:
%   cfg = ARC_makeDecodingCfg('C:\Work\ARC\ARC');

% ── basic paths ──────────────────────────────────────────────────────────
cfg.root         = rootDir;
cfg.maskFile     = 'ARC3_anatgw.nii';
cfg.fMaskFile    = 'ARC3_fanatgw3_pos.nii';

% ── analysis parameters ─────────────────────────────────────────────────
cfg.bins         = 7;              % behavioural bins
cfg.binCentre    = (cfg.bins+1)/2; % middle bin index (4 when bins = 7)
cfg.subjects     = [1 2 3];        % subject IDs to analyse
cfg.behavVector  = 2;              % index into behav.ratings (:,v_id)

cfg.popChoice    = 'pos';          % 'pos','neg','mut'
cfg.doZscore     = false;          % voxel-wise z-score?
cfg.noisePoolCut = false;          % drop GLMdenoise noisepool?
cfg.sigVoxCut    = false;          % keep only sig. odor voxels?
cfg.shuffleBehav = false;          % sanity shuffle

% ── ROIs ────────────────────────────────────────────────────────────────
% cfg.ROIs.names = {'Insula','Hipp','DLPFC','A1','wm'};
% cfg.ROIs.masks = {'rwinsula.nii','rwHipp.nii','rwDLPFC.nii',...
                  % 'rwAud.nii','rwm_main.nii'};
cfg.ROIs.names = {'PC','AMY','OFC','VMPFC'};
cfg.ROIs.masks = {'rwPC.nii','rwAmygdala.nii','rwofc.nii','rwvmpfc.nii'};

% ── save / load ─────────────────────────────────────────────────────────
cfg.saveDir      = fullfile(rootDir,'Decoding','RSA_correct','Basic_extraROI');
cfg.RSAfile      = fullfile(rootDir,'results','temp_main_3','ARC_wts.mat'); % w_score_mat etc.
cfg.behavFile    = fullfile(rootDir,'supporting_files','NEMO_perceptual2.mat');
cfg.trialFile    = fullfile(rootDir,'supporting_files','valence_trials.mat');
cfg.libsvmPath   = fullfile(rootDir,'supporting_files','libsvm');  % if needed

cfg.nfold = 4;
% ── helper derived values ───────────────────────────────────────────────
cfg.nbSubjects   = numel(cfg.subjects);
cfg.nbROIs       = numel(cfg.ROIs.names);
end
