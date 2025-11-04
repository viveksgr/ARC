function cfg = ARC_makeDecodingCfg(rootDir,seed)
%ARC_makeDecodingCfg  Create a struct with all analysis settings
%
% Usage:
%   cfg = ARC_makeDecodingCfg('C:\Work\ARC\ARC');

if nargin<2
    seed = [];
end

% Freq toggles
cfg.samepop = false;
cfg.popChoice    = 'mut';          % 'pos','neg','mut'
% ── save / load ─────────────────────────────────────────────────────────
cfg.saveDir      = fullfile(rootDir,'Decoding','searchl_test');
cfg.group_model = false;
cfg.searchl = true;
% cfg.splithalf_ = false;
cfg.rget = 2;
cfg.disc =false;
% cfg.ig_neut = false;
cfg.includeneut = true;

cfg.sesswise = true;
% cfg.RSAfile      = fullfile(rootDir,'results\mixedmaster100','wt_mats.mat'); % w_score_mat
cfg.nfold = 4;
cfg.prctile = true; 
cfg.thrs_1 = 0.7;
cfg.thrs_2 = 0.5;
% cfg.thrs_1 = 0.3;
% cfg.thrs_2 = 0.1;
cfg.cross = true;
cfg.mapper = true;


% ── basic paths ──────────────────────────────────────────────────────────
cfg.root         = rootDir;
cfg.maskFile     = 'ARC3_anatgw.nii';
cfg.fMaskFile    = 'ARC3_fanatgw3_pos.nii';

% ── analysis parameters ─────────────────────────────────────────────────
cfg.bins         = 7;              % behavioural bins
cfg.binCentre    = (cfg.bins+1)/2; % middle bin index (4 when bins = 7)
cfg.subjects     = [1 2 3];        % subject IDs to analyse
cfg.behavVector  = 2;              % index into behav.ratings (:,v_id)

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
cfg.histfile   = fullfile(rootDir,'supporting_files','histsplit_ids_mat.mat');
cfg.seed = seed;

% cfg.RSAfile      = fullfile(rootDir,'results','temp_main_2','ARC_RSA.mat'); % w_score_mat etc.
cfg.behavFile    = fullfile(rootDir,'supporting_files','NEMO_perceptual2.mat');
cfg.trialFile    = fullfile(rootDir,'supporting_files','valence_trials.mat');
cfg.libsvmPath   = fullfile(rootDir,'supporting_files','libsvm');  % if needed

% ── helper derived values ───────────────────────────────────────────────
cfg.nbSubjects   = numel(cfg.subjects);
cfg.nbROIs       = numel(cfg.ROIs.names);
cfg.splithalf_  = true;
cfg.nested = false;
cfg.verbose = true;

end
