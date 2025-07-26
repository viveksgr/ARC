function res = ARC_runDecoding(cfg)
%ARC_runDecoding  Main entry point for the decoding analysis

% addpath(cfg.libsvmPath);                              % libsvm, if required
% behav = load(cfg.behavFile);                          % behavioural ratings
trfs = load(cfg.trialFile);

% containers for results
beta   = zeros(cfg.nbSubjects,cfg.nbROIs,2);
tstat  = beta;
dfmat  = beta;

for sIdx = 1:cfg.nbSubjects
    subj = cfg.subjects(sIdx);
    fprintf('Subject %02d\n',subj);

    % ▸ behavioural vector for this subject
    % beh = behav.behav(subj).ratings(:,cfg.behavVector);
    beh = trfs.behav_ratings_{sIdx};
    if cfg.shuffleBehav, beh = beh(randperm(numel(beh))); end
    beh = normalize(beh,'medianiqr');

    % ▸ per-ROI loop
    for r = 1:cfg.nbROIs
        [neural,posIdx,negIdx,mutIdx,labels] = ...
            ARC_prepareROI(cfg,subj,r,beh);

        if cfg.ig_neut
            labels_pos = labels>cfg.binCentre;
            labels_neg = labels<cfg.binCentre;
        else
            labels_neg = selectLess4AndHalf4(labels);
            labels_pos = ~labels_neg;
        end

        if cfg.splithalf_
            group_vec = trfs.gr{subj};
            trial_id_mat = load(cfg.histfile);
            cfg.odor_select = trial_id_mat.idxmat(:,cfg.seed);
            trial_inc = ~ismember(group_vec,cfg.odor_select); % Choose trials NOT used in RSA
            labels_pos = and(labels_pos,trial_inc);
            labels_neg = and(labels_neg,trial_inc);
        end

        % choose voxel population ---------------------------------------
        if cfg.samepop

            switch cfg.popChoice
                case 'pos'
                    X1 = neural(labels_pos,posIdx);   % V+ coded by positive pop
                    X2 = neural(labels_neg,posIdx);   % V- coded by negative pop
                    y1 = labels(labels_pos);
                    y2 = labels(labels_neg);
                case 'neg'
                    X1 = neural(labels_pos,negIdx);   % V+ coded by negative pop
                    X2 = neural(labels_neg,negIdx);   % V- coded by positive pop
                    y1 = labels(labels_pos);
                    y2 = labels(labels_neg);
                case 'mut'                                          % mutually tuned voxels
                    X1 = neural(labels_pos,mutIdx);
                    X2 = neural(labels_neg,mutIdx);
                    y1 = labels(labels_pos);
                    y2 = labels(labels_neg);
                otherwise, error('Bad popChoice');
            end

            if cfg.nested
                % ▸ decode two tasks (valence+, valence-) -----------------------
                [r1,~,t1] = ARC_regress_nested2_normed(X1,y1,cfg.nfold,1);
                [r2,~,t2] = ARC_regress_nested2_normed(X2,y2,cfg.nfold,1);
            else
                % ▸ decode two tasks (valence+, valence-) -----------------------
                [r1,~,t1] = ARC_regress_normed(X1,y1,cfg.nfold,1);
                [r2,~,t2] = ARC_regress_normed(X2,y2,cfg.nfold,1);
            end


            % store
            beta  (sIdx,r,:) = [r1,r2];
            tstat (sIdx,r,:) = [t1,t2];
            dfmat (sIdx,r,:) = [length(y1),length(y2)];

        else


            switch cfg.popChoice
                case 'pos'
                    X1p = neural(labels_pos,posIdx);   % V+ coded by positive pop
                    X2p = neural(labels_neg,posIdx);   % V- coded by positive pop
                    X1n = neural(labels_pos,negIdx);   % V+ coded by negative pop
                    X2n = neural(labels_neg,negIdx);   % V- coded by negative pop
                    y1 = labels(labels_pos);
                    y2 = labels(labels_neg);


                    [ r1,~,t1] = ARC_regress_wrapper(X1p,y1,X2p,y2, cfg.nfold,1,cfg.nested);
                    [ r2,~,t2] = ARC_regress_wrapper(X2p,y2,X1p,y1, cfg.nfold,1,cfg.nested);

                case 'neg'
                    X1p = neural(labels_pos,posIdx);   % V+ coded by positive pop
                    X2p = neural(labels_neg,posIdx);   % V- coded by positive pop
                    X1n = neural(labels_pos,negIdx);   % V+ coded by negative pop
                    X2n = neural(labels_neg,negIdx);   % V- coded by negative
                    y1 = labels(labels_pos);
                    y2 = labels(labels_neg);

             [  r1,~,t1] = ARC_regress_wrapper(X1n,y1,X2n,y2, cfg.nfold,1,cfg.nested);
             [ r2,~,t2] = ARC_regress_wrapper(X2n,y2,X1n,y1, cfg.nfold,1,cfg.nested);
             

                case 'mut'                                          % mutually tuned voxels
                    X1 = neural(labels_pos,mutIdx);
                    X2 = neural(labels_neg,mutIdx);
                    y1 = labels(labels_pos);
                    y2 = labels(labels_neg);

      
             [  r1,~,t1] = ARC_regress_wrapper(X1,y1,X2,y2, cfg.nfold,1,cfg.nested);
             [ r2,~,t2] = ARC_regress_wrapper(X2,y2,X1,y1, cfg.nfold,1,cfg.nested);
             


                otherwise, error('Bad popChoice');
            end
     
            % store
            beta  (sIdx,r,:) = [r1,r2];
            tstat (sIdx,r,:) = [t1,t2];
            dfmat (sIdx,r,:) = [length(y1),length(y2)];

        end
    end
end

% ▸ group-level statistics & plotting -----------------------------------
% ARC_plotDecoding(beta,tstat,dfmat,cfg);
res.beta = beta;
res.dfmat = dfmat;


if cfg.verbose
    mkdir(cfg.saveDir);
    [pvals] = ARC_decoding_pvals(beta,dfmat);     % your existing helper
    ARC_barplot_sig(beta,pvals)
    xticks(1:cfg.nbROIs)
    xticklabels(cfg.ROIs.names);
    savefig(fullfile(cfg.saveDir,"ARC_decoding.fig"))
    print(fullfile(cfg.saveDir,'ARC_decoding'),'-dpng')
    print(gcf,'-vector','-dsvg',[fullfile(cfg.saveDir,'ARC_decoding'),'.svg']) % svg
    save(fullfile(cfg.saveDir,'ARC_decoding.mat'),'beta','tstat','dfmat','cfg');
end
end
