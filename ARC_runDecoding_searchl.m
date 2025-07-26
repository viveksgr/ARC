function res = ARC_runDecoding_searchl(cfg)
%ARC_runDecoding  Main entry point for the decoding analysis

% addpath(cfg.libsvmPath);                              % libsvm, if required
% behav = load(cfg.behavFile);                          % behavioural ratings
tic
trfs = load(cfg.trialFile);

% containers for results
beta   = cell(cfg.nbSubjects,cfg.nbROIs,2);
dfmat  = zeros(cfg.nbSubjects,cfg.nbROIs,2);
mkdir(cfg.saveDir);
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
        fprintf('Anat ROI %02d\n',r);

        [neural,~,~,~,labels] = ...
            ARC_prepareROI(cfg,subj,r,beh);
        assert(sum(isnan(neural),'all')==0,'Nan elements in neural matrix')
        anatmask = ARC_configROI_mask(cfg,subj,r);
        assert(sum( anatmask,'all')==size(neural,2),'Mismatch between mask and ROI')


        [dec_val1,dec_val2,df] = ARC_searchlightdec_ROI(cfg,anatmask,neural,labels,subj,r);
        beta{subj,r,1} = dec_val1;
        beta{subj,r,2} = dec_val2;
        dfmat(subj,r,1) = df(1);
        dfmat(subj,r,2) = df(2);

    end
end

% ▸ group-level statistics & plotting -----------------------------------
% ARC_plotDecoding(beta,tstat,dfmat,cfg);
res.beta = beta;
res.dfmat = dfmat;
save(fullfile(cfg.saveDir,'ARC_decoding.mat'),'res','cfg');
toc

[beta_mat,pvals] = ARC_searchlight_transform(beta,dfmat);
pvals = pvals+0.1;
ARC_barplot_sig(beta_mat,pvals)
xticks(1:cfg.nbROIs); xticklabels(cfg.ROIs.names);
ylabel('% sig searchlights'); legend({'Valence','Salience'});
ylim([0 25])
% title(sprintf('Decoding: %s population',cfg.popChoice));
print(fullfile(cfg.saveDir,'ARC_decoding'),'-dpng');
savefig(fullfile(cfg.saveDir,'ARC_decoding'));
% ARC_barplot_sig(beta_mat,pvals)

%% Combine NII across ROIs
thrs = zeros(3,4);
for ss = 1:3
    for  k_id = [1 2]

        fname = sprintf('ARC%02d_val%01d',ss,k_id);
        flist = dir(fullfile(cfg.saveDir,'ROIwise',sprintf('%s*.nii',fname)));
        ffile = {}; for ig = 1:numel(flist); ffile{ig} = fullfile(flist(ig).folder, flist(ig).name); end
        fdata = []; for ig = 1:numel(flist); fdata(:,:,:,ig) = spm_read_vols(spm_vol(ffile{ig})); end
        fdata(fdata==0) = nan;
        fmean = mean(fdata,4,'omitmissing');

        rdata = fmean(~isnan(fmean));
        df = dfmat(ss,1,k_id);

        [~,  thr_pos, ~, thr_neg] = ARC_bidirectionalFDR(rdata, df);
        kid_in = 2*(k_id-1);
        thrs(ss,  kid_in +1) =  r2t(thr_pos,df);
        thrs(ss,  kid_in +2) =  r2t(thr_neg,df);

        subjDir  = fullfile(cfg.root,sprintf('ARC%02d',ss),'single');

        write_reshaped_nifty(fmean, cfg.saveDir, false, fullfile(subjDir,cfg.maskFile), sprintf('ARC%02d_val%01d',ss,k_id ));
        write_reshaped_nifty(-fmean, cfg.saveDir, false, fullfile(subjDir,cfg.maskFile), sprintf('nARC%02d_val%01d',ss,k_id ));

    end
end
res.thrs = thrs;

%% Normalize and add
% k_id = [val1 sal1 val2 sal2];

k_thrs = zeros(1,4);
for k_id = [1,2];
    kid_in = 2*(k_id-1);
    fname = sprintf('wARC*_val%01d',k_id);
    flist = dir(fullfile(cfg.saveDir,sprintf('%s*.nii',fname)));
    ffile = {}; for ig = 1:numel(flist); ffile{ig} = fullfile(flist(ig).folder, flist(ig).name); end
    rmaps= []; for ig = 1:numel(flist); rmaps(:,:,:,ig) = spm_read_vols(spm_vol(ffile{ig})); end
    Ns = squeeze(dfmat(:,1,k_id));
    [meanR,k_thrs(kid_in +1) ,k_thrs(kid_in +2)] = combineRmaps_meta_bi_weighted(rmaps, Ns);
    write_reshaped_nifty(meanR, cfg.saveDir, true, fullfile( cfg.saveDir,'wARC01_val1.nii'), sprintf('ARCmain_val%01d',k_id ));
    write_reshaped_nifty(-meanR, cfg.saveDir, true, fullfile( cfg.saveDir,'wARC01_val1.nii'), sprintf('nARCmain_val%01d',k_id ));
end
res.metathr = k_thrs;

%% Bar plots
thr_mat = zeros(3,4,2);
for ii = 1:3
    for jj = 1:2
        thr_mat(ii,:,jj) = thrs(ii,jj)*ones(4,1);
    end
end

[beta_mat] =  ARC_thresholdPercentages(beta(:,:,1),thr_mat);
pvals = pvals+0.1;
ARC_barplot_sig(beta_mat,pvals)
xticks(1:cfg.nbROIs); xticklabels(cfg.ROIs.names);
ylabel('% sig searchlights'); legend({'Valence','Salience'});
% title(sprintf('Decoding: %s population',cfg.popChoice));
print(fullfile(cfg.saveDir,'ARC_decoding'),'-dpng');
savefig(fullfile(cfg.saveDir,'ARC_decoding'));
ARC_barplot_sig(beta_mat,pvals)


%     % Image avg.
% matlabbatch{1}.spm.util.imcalc.input = {
%                                         'C:\Data\ARC\ARC02\mediation\Searchlightx1_4anat_median\wRSAx1OFC_seedPC_valpos.nii,1'
%                                         'C:\Data\ARC\ARC01\mediation\Searchlightx1_4anat_median\wRSAx1OFC_seedPC_valpos.nii,1'
%                                         'C:\Data\ARC\ARC03\mediation\Searchlightx1_4anat_median\wRSAx1OFC_seedPC_valpos.nii,1'
%                                         };
% matlabbatch{1}.spm.util.imcalc.output = 'OFC_PC_val';
% matlabbatch{1}.spm.util.imcalc.outdir = {''};
% matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3)/3';
% matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
% matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
% matlabbatch{1}.spm.util.imcalc.options.mask = 0;
% matlabbatch{1}.spm.util.imcalc.options.interp = 1;
% matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
% spm_jobman('run', matlabbatch);
%

% if cfg.verbose
%     mkdir(cfg.saveDir);
%     [pvals] = ARC_decoding_pvals(beta,dfmat);     % your existing helper
%     ARC_barplot_sig(beta,pvals)
%     xticks(1:cfg.nbROIs)
%     xticklabels(cfg.ROIs.names);
%     savefig(fullfile(cfg.saveDir,"ARC_decoding.fig"))
%     print(fullfile(cfg.saveDir,'ARC_decoding'),'-dpng')
%     print(gcf,'-vector','-dsvg',[fullfile(cfg.saveDir,'ARC_decoding'),'.svg']) % svg
%     save(fullfile(cfg.saveDir,'ARC_decoding.mat'),'beta','tstat','dfmat','cfg');
% end
end
