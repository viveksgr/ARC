function results = ARC_run_RSA_full(cfg)
% RUN_RSA_FULL  loops over subjects & ROIs and collects everything
%
%   results has fields:
%       .wSub       {S x R}   beta weights
%       .tSub       {S x R}   t-scores
%       .pSub       [S x R x 2]  permutation p
%       .groupStats            (struct with FDR / bar‐plot data)

nSub   = numel(cfg.subjectList);
nROI   = numel(cfg.anatNames);

% containers
wSub = cell(nSub,nROI);
tSub = cell(nSub,nROI);
pSub = cell(nSub,nROI);
pcorr = zeros(nSub,nROI);
w_score_mat = cell(nSub,nROI);
w_score_mat_dist = cell(nSub,nROI);
w_mat = cell(nROI,1);

for si = 1:nSub
    sid = cfg.subjectList(si);
    fprintf('\nSubject %02d\n',sid);

    subjDat = ARC_load_subject_data(sid,cfg);                 % behav + GLMs
    for ri = 1:nROI
        roiName = cfg.anatNames{ri};
        fprintf('  ROI: %s\n',roiName);

        roiDat = ARC_extract_roi_trials(subjDat,ri,cfg,si);      % mask+voxel list
       
          % if (ri==3)
          %       'error'
          % end

        [w,wt_dist,res] = ARC_run_RSA_for_subject(roiDat,cfg);        % <-- core statistics
        wSub{si,ri}     = w;
        tSub{si,ri}     = wt_dist;
        pSub{si,ri}     =  res;
        
        if ~cfg.runpcmcollapse
        pcorr(si,ri) = res.t_corr;
        % w_score_mat{si,ri} = res.wt_mat;
        w_score_mat_dist{si,ri} = res.w_scores;
        w_mat{ri} = cat(1,w_mat{ri},res.w_scores);
        end
    end
end

% Save results
results.wSub = wSub;
results.tSub = tSub;
results.pSub = pSub;
results.w_score_matdist = w_score_mat_dist;

% D0 stats and plots
% ---------- combine & plot ----------
if cfg.verbose
group = ARC_make_group_stats(wSub, tSub,cfg);
results.group = group;
create_ARCscatter(w_score_mat_dist,cfg)
end

save(fullfile(cfg.saveRoot,'ARC_RSA.mat'));
fprintf("Done")
end
