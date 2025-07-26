function cfg = ARC_make_default_config()
% Toggle between fullodor set and partial: toggle; name/dir, nshuffle,
% verbose and  split_half


% ---------- master toggles ----------
cfg.valenceSplit = false;           % val+/val–
cfg.splithalf_   = false;
cfg.voxwiseshuff = true;
cfg.usepcm = false;
cfg.runpcmcollapse = false;
cfg.runvoxwise   = false;
cfg.plsregress = true;

cfg.splitneupop = 'all'; % 'neg';'neut' or 'all';
cfg.numCtrl      = true;
cfg.sizeCtrl     = false;
cfg.intensityReg = false;
cfg.sniffCtrl    = false;
cfg.zscoreRows   = true;
cfg.runShuffle   = false;
cfg.runSqDist    = false;
cfg.shufftest = false;
cfg.numCtrlbin = 200;

cfg.ignore_neut = false;
cfg.plotlims = [-0.2 1];
% cfg.seed = 1; % put cfg seed as 1-100

% ---------- analysis basics ----------
cfg.mainRoot   = 'C:\Work\ARC\ARC';
cfg.modelName  = 'LDA_temp_shuff';
cfg.saveRoot   = fullfile(cfg.mainRoot,'results',cfg.modelName);

cfg.subjectList = 1:3;            % which subject IDs to analyse
cfg.anatNames   = {'PC','AMY','OFC','VMPFC'};
cfg.anatMasks   = {'rwPC.nii','rwAmygdala.nii','rwofc.nii','rwvmpfc.nii'};
% cfg.anatNames = {'Insula','Hipp','DLPFC','A1','wm'};
% cfg.anatMasks = {'rwinsula.nii','rwHipp.nii','rwDLPFC.nii','rwAud.nii','rwm_main.nii'};

cfg.behavFile = fullfile(cfg.mainRoot,'supporting_files','NEMO_perceptual.mat');
cfg.maskFile  = 'ARC3_anatgw.nii';
cfg.fMaskFile = 'ARC3_fanatgw3_pos.nii';

cfg.nBin      = 7;
cfg.nShuffle  = 1000;
cfg.binzPart1 = ceil((cfg.nBin+1)/2);   % helper for val+/val–
cfg.binzPart2 = cfg.nBin - cfg.binzPart1 + 1;
cfg.v_ids = [2 2 2];
cfg.nodor = 160;

cfg.histfile = fullfile(cfg.mainRoot,'supporting_files','histsplit_ids_mat.mat');
cfg.verbose = true; % save results figs etc
cfg.thrs_1 = 0.3;
cfg.thrs_2 = 0.1;

% % Confirmation
% if cfg.voxwiseshuff
%     fprintf('Voxwise shuffling control analysis')
% end
