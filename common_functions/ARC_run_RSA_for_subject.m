 function [w,wt_dist,res] = ARC_run_RSA_for_subject(roiDat,cfg)

% Build behavioural RSM(s)
[val1, val2, maskMat,v1,v2] = ARC_build_behav_RSM(cfg);

% Bin neural RSM
behav_ratings = roiDat.behav.ratings(:,cfg.v_ids(roiDat.si));
behav_ratings = normalize(behav_ratings,'medianiqr');

if cfg.splithalf_
    group_vec = roiDat.group_vec;
    trial_id_mat = load(cfg.histfile);
    cfg.odor_select = trial_id_mat.idxmat(:,cfg.seed);
    trial_inc = ismember(roiDat.group_vec,cfg.odor_select);
  else
    trial_inc = ismember(roiDat.group_vec,roiDat.group_vec);
end

behav_ratings_exp = behav_ratings(roiDat.group_vec(trial_inc));

M = roiDat.betatrials';
idx_anat = ARC_classifyVoxelSign(M,0.05,'t');
if strcmp(cfg.splitneupop,'pos')
    sel = idx_anat==1;
    fprintf('Num vox: %02d\n',sum(sel))
elseif strcmp(cfg.splitneupop,'neg')
    sel = idx_anat==-1;
    fprintf('Num vox: %02d\n',sum(sel))
elseif strcmp(cfg.splitneupop,'neut')
    sel = idx_anat==0;
    fprintf('Num vox: %02d\n',sum(sel))
elseif strcmp(cfg.splitneupop,'all')
    sel = idx_anat<99;
else
    error('Wrong pop description')
end

[modelmd_binned,modelmd_binned_mat,bin] = ARC_binAndTransform_numctrl(roiDat.betatrials(sel,trial_inc), behav_ratings_exp, cfg.nBin, cfg.numCtrlbin,cfg.nShuffle);

modelmd_binned_mat = permute(modelmd_binned_mat,[3 1 2]);
[~,modelmd_binned_shuff] = ARC_binAndTransform_shuffcoarse(roiDat.betatrials(sel,trial_inc), behav_ratings_exp, cfg.nBin, cfg.numCtrlbin,cfg.nShuffle);
modelmd_binned_shuff = permute(modelmd_binned_shuff,[3 1 2]);

if cfg.shufftest
    modelmd_binned_mat = modelmd_binned_shuff;
end

% % Fit multiple-regression RSA
% [w,t] = ARC_multicomputeWeights_tsc([val1(maskMat) val2(maskMat)],  modelmd_corrcoef(maskMat));

wt_dist = zeros(cfg.nShuffle,3);
for zz = 1:cfg.nShuffle % shuffle
    modelmd_binneds = squeeze(modelmd_binned_mat(zz,:,:));
    if  cfg.zscoreRows
        modelmd_binneds = zscore(modelmd_binneds,[],2);
    end
    modelmd_corrcoef = corrcoef(modelmd_binneds);  
    [wt_dist(zz,:), ~] = ARC_multicomputeWeights_tsc([val1(maskMat) val2(maskMat)],modelmd_corrcoef(maskMat));
end
w = mean(wt_dist)';

% Shuffle test
wt_dist = zeros(cfg.nShuffle,3);
for zz = 1:cfg.nShuffle % shuffle
    modelmd_binneds = squeeze(modelmd_binned_shuff(zz,:,:));
    if  cfg.zscoreRows
        modelmd_binneds = zscore(modelmd_binneds,[],2);
    end
    modelmd_corrcoef = corrcoef(modelmd_binneds);  
    [wt_dist(zz,:), ~] = ARC_multicomputeWeights_tsc([val1(maskMat) val2(maskMat)],modelmd_corrcoef(maskMat));
end

% Voxelwise analysis
if cfg.runvoxwise
    res = ARC_voxwise(val1, val2, modelmd_binned_mat,cfg);
    % M = roiDat.betatrials(:,trial_inc);
    % runIdx = blockRowIndices(roiDat.trialstruct);   
    % 
    % % tic
    % if cfg.voxwiseshuff
    % M = M(:,randperm(size(M,2)));
    % end
    % 
    % % if cfg.valenceSplit
    %     vn = v1;
    %     vp = v1;
    %     vn(vn>0)=0;
    %     vp(vp<0)=0;
    %     res = ARC_pcm(vn,vp,M,bin,runIdx(trial_inc),cfg);
    % else
    %     res = ARC_pcm(v1,v2,M,bin,runIdx(trial_inc),cfg);   
    % end
    % if and(cfg.usepcm,~cfg.runpcmcollapse)
    %     w(2:3) = res.beta;
    % end
    %  res3 = ARC_pcm(vn,vp,M_shuff,bin,runIdx(trial_inc));
    % toc
else
   res = ARC_voxwise(val1, val2, modelmd_binned_mat(1:2,:,:),cfg);
end
res.df = sum(sel);

if cfg.voxwiseshuff
    res_shuff = ARC_voxwise(val1, val2, modelmd_binned_shuff ,cfg);
    res.res_shuff = res_shuff;
end

if cfg.plsregress
    X = roiDat.betatrials(sel,trial_inc)';
    if cfg.voxwiseshuff
        X = X(randperm(size(X,1)),:);
    end

    % vn = v1;
    % vp = v1;
    % vn(vn>0)=0;
    % vp(vp<0)=0;
    % v1_vec = vp(bin)';
    % v2_vec = vn(bin)';
    % res = ARC_voxelPLS(X, v1_vec, v2_vec);

    vn = v1;
    vp = v1;
    vn(vn>0)=0;
    vp(vp<=0)=0;
    vn = logical(vn);
    vp = logical(vp);
    v1_vec = vp(bin)';
    v2_vec = vn(bin)';
    tic
    try
    [res] = ARC_twoAxisLDA(X, v1_vec, v2_vec);
    catch
        'beep'
    end
    toc
    % vn = v1;
    % vp = v1;
    % vn(vn>0)=0;
    % vp(vp<0)=0;
    % vn = ~logical(vn);
    % vp = ~logical(vp);
    % v1_vec = vn(bin)';
    % v2_vec = vp(bin)';
    % 
    % v11 = logical([0 0 0 0 0 0 1 1]);
    % v22 = logical([1 1 0 0 0 0 0 0]);
    % [res] = ARC_fourAxisLDA(X, v1_vec, v2_vec,v11(bin),v22(bin));
end

% if ~cfg.verbose
%     res.wt_mat = [];
% end
end
