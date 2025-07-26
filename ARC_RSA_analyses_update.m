% %% RSA Analyses
% This script performs Representational Similarity Analysis (RSA) on fMRI data
% to investigate valence, salience, appetitive and aversive representations in neural activity.

% Key Features:
% - Loads behavioral and fMRI data for multiple subjects.
% - Constructs representational similarity matrices (RSMs) for valence and salience.
% - Correlates behavioral RSMs with neural response RSMs.
% - Implements permutation testing for statistical significance of RSA weights.
% - Computes voxel-wise and region-level statistics, with FDR correction.
% - Visualizes results with bar plots, similarity matrices, and scatter density plots.
% - Saves RSA weights and voxel-wise maps as figures and NIfTI files.

% Inputs:
% - Behavioral data from NEMO_perceptual2.mat (pleasantness ratings of odors).
% - Anatomical and functional masks for ROI-based analysis.
% - Single trial neural responses from GLM single package.

% Outputs:
% - RSA statistics: Weights, t-scores, permutation p-values.
% - Figures: RSA bar plots, voxel proportions, scatter density plots.
% - NIfTI files: Voxel-wise RSA maps for positive and negative valence.

% Usage:
% Adjust analysis settings such as `binz`, `anat_names`, and `dirs` for your data.
% Run the script to compute RSA, visualize results, and save outputs.

% VivekSagar2016@u.northwestern.edu;
% January 2025

%% General Settings
% Main settings adjusted:
% val_sep flips whether the RSA is for valence/salience (false) or
% appetitive/aversive domain (true)

% The script maybe in demo-mode (demomode = true) that loads trial-averaged data instead of
% single trial data. Adjust it to false for the full model training.

demomode = false; % Current demo is only for valsep, num_cntrl, sz_cntrl, intens_reg = false;
valsep = false;
num_cntrl = true; % Control analysis for bin size
numctrl_bin = 100;

sz_ctrl = false; % Control analysis for ROI size
intens_reg = false; % Control analysis for Intensity
sniff_ctrl = false;
bias_correct = false;
seq = false; % Regressout valp valn
run_shuff = false;
run_sqdist = false;
run_voxwise = false;

% root = 'C:\Work\ARC\Scripts';
mainroot = 'C:\Work\ARC\ARC';
modelname ='main_oldscript_100bin';
savepath = fullfile(mainroot,'results',modelname) ;
% v_ids = [1 1 1];
v_ids = [2 2 2];
% v_ids = [14 13 13]; % Warm
% v_ids = [0 18 18];

maskfile =  'ARC3_anatgw.nii';
fmaskfile = 'ARC3_fanatgw3_pos.nii';
fmasker = false;
binz =7 ; % Number of bins
if mod(binz,2)==0; binzpart1 = binz/2; binzpart2 = binzpart1+1; else; binzpart1 = (binz+1)/2 ; binzpart2 = binzpart1; end

% anat_names = {'PC','AMY','OFC','VMPFC','Insula','Hipp','DLPFC','A1','wm'};
% anat_masks = {'rwPC.nii','rwAmygdala.nii','rwofc.nii','rwvmpfc.nii','rwinsula.nii','rwHipp.nii','rwDLPFC.nii','rwAud.nii','rwm_main.nii'};

anat_names = {'PC','AMY','OFC','VMPFC'};
anat_masks = {'rwPC.nii','rwAmygdala.nii','rwofc.nii','rwvmpfc.nii'};

% anat_names = {'Insula','Hipp','DLPFC','A1','wm'};
% anat_masks = {'rwinsula.nii','rwHipp.nii','rwDLPFC.nii','rwAud.nii','rwm_main.nii'};
nanat = length(anat_names);

single_n = false; % Noisepool from GLM single
single_c = false;  % Voxelwise cutoff
zscorer = true; % Normalization
noblock = false;

nodor = 160;
nshuff = 1000;
behav = load(fullfile(mainroot,'supporting_files','NEMO_perceptual.mat'));
if sniff_ctrl; load(fullfile(mainroot,'supporting_files','sniff_corr_ctrl_correct.mat')); end

% load(fullfile(statpath,'fir_cv.mat'))
fprintf('\n')

% Result matrics
rsa_P1 = nan(3,nanat,2); % For t-scores in anatomical areas
rsa_P1wt = nan(3,nanat,2); % For beta-wts in anatomical areas
rsa_Pp = nan(3,nanat,2); % p-values of RSA
rsa_Pps = nan(3,nanat,nshuff,2); % Null distributions
rsa_Pcorr = nan(3,nanat); % Correlation of beta weights across domains
rsa_prop = nan(3,nanat,3);  % Proportion of significant voxels in each regions
modelbinned_mat = cell(3,nanat); % Matrix of subjectxROI binned values of neural responses for demo

t_score_mat = cell(3,4);
w_score_mat = cell(3,4);
w_mat = cell(nanat,1);
% FDR corrected Stats
thr_fdr = zeros(3,2);
thr_fdranat = zeros(3,2,nanat);

% Figure across subjects
figure1 = figure('OuterPosition',[297 183 1209 737]);
hold on
kk = 1;
for s = [1 2 3] % Subject
    v_id = v_ids(s);
    fprintf('Subject: %02d\n',s)

    if v_id==0
        continue
    else

        anatdir = fullfile(mainroot,sprintf('ARC%02d',s),'single');
        if demomode;  anatdir = fullfile(mainroot,'supporting_files',sprintf('ARC%02d',s)); load(fullfile(mainroot,'supporting_files','modelbinned_mat.mat')); end

        % Construction of median split
        behav_ratings = behav.behav(s).ratings(:,v_id);

        % fprintf('Shuffled ratings in progress')
        % behav_ratings = behav_ratings(randperm(length(behav_ratings)));
        % behav_ratings = behav_ratings-median(behav_ratings);

        behav_int = behav.behav(s).ratings(:,1);
        if ~bias_correct; behav_ratings = normalize(behav_ratings,'medianiqr'); end
        ms1 = min(behav_ratings);
        ms2 = max(behav_ratings);
        behav_int = normalize( behav_int,'medianiqr');
        md = 0;

        % Gray Matter, Functional and Anatomical Masks
        mask = (spm_read_vols(spm_vol(fullfile(anatdir, maskfile)))); % Mask used to construct odor files
        mask(isnan(mask))=0;
        mask = logical(mask);
        fmask = (spm_read_vols(spm_vol(fullfile(anatdir, fmaskfile)))); % Mask used to examine voxels in RSA
        fmask(isnan(fmask))=0;
        if ~fmasker
            fmask = fmask +0.1;
        end
        fmask = logical(fmask); % Only choose voxels with significant odor evoked activity
        fmask_1d = fmask(mask);
        marea = and(fmask,mask);

        % Model names
        masks_set = [];
        masks_set_cell = {};
        anatmasks = [];
        for ii = 1:length(anat_masks)
            m1 = spm_read_vols(spm_vol(fullfile(anatdir,anat_masks{ii})));
            m1(isnan(m1))=0;
            m1(m1<=0.01)=0;
            m1(m1>0) = 1;
            m1 = m1(mask);
            masks_set(:,ii)=m1(fmask_1d);
            fprintf('area count: %04d\n',sum(m1(fmask_1d)))

            masks_set_cell{ii} = fmask_1d(logical(m1));
            anatmask = (spm_read_vols(spm_vol(fullfile(anatdir, anat_masks{ii}))));
            anatmask(isnan(anatmask)) = 0;
            anatmask = logical(anatmask);
            anatmask = and(anatmask,marea);
            anatmasks(:,:,:,ii) = anatmask;
        end
        anat_cell{s}= anatmasks;
        masks_set(isnan(masks_set))=0;

        linux_config = false;
        warning('off','all')
        map_area = zeros([size(anat_cell{s}) 2]);

        %% Rep similarity
        S_mat = zeros(length(anat_names),2);
        S_mat2 = zeros(length(anat_names),2);

        onsets = load(fullfile(anatdir,sprintf('conditions_NEMO%02d.mat',s)),'onsets');
        onsets = onsets.onsets;
        group_vec = cell(nodor,1);
        unity = [];
        for ii2 = 1:nodor
            group_vec{ii2} = ii2*ones(length(onsets{ii2}),1);
            unity = blkdiag(unity,ones(length(onsets{ii2})));
        end
        group_vec = vertcat(group_vec{:});
        [~,argsort] = sort(vertcat(onsets{:}));
        group_vec = group_vec(argsort);
        unity = unity(argsort,argsort);
        utl_mask = logical(triu(ones(length(unity)),1)); % All possible odors

        behav_ratings_ = behav_ratings(group_vec);
        if bias_correct
            bias_idx = ARC_median_balance_zero( behav_ratings_);
            behav_ratings_ = behav_ratings_(bias_idx);
        end

        for ii = 1:length(anat_names)

            fprintf('area:%02d\n',ii)

            if ~demomode % Extract raw single trial responses (not functional in the demo version)
                modelmd_ = load(fullfile(anatdir,anat_names{ii},'TYPED_FITHRF_GLMDENOISE_RR.mat'),'modelmd','noisepool');
                % modelmd_ = load(fullfile(anatdir,anat_names{ii},'TYPEB_FITHRF.mat'),'modelmd','noisepool');

                modelmd = squeeze(modelmd_.modelmd);

                if single_c
                    noisepool = modelmd_.noisepool;
                    modelmd = modelmd(masks_set_cell{ii},:);
                    noisepool = noisepool(masks_set_cell{ii});
                end
                if single_n
                    noisepool = modelmd_.noisepool;
                    modelmd2 = modelmd(~noisepool,:);
                else
                    modelmd2 = modelmd;
                end
                fprintf('size:%02d\n',size(modelmd,1))
                [r1,~] = find(isnan(modelmd2));
                modelmd2(r1,:) = [];
            end

            if bias_correct
                modelmd2 = modelmd2(:, bias_idx);
            end

            if ~num_cntrl
                % Basic model
                if ~demomode
                    modelmd_binned = ARC_binAndTransform(modelmd2, behav_ratings_, binz, [ms1 ms2]);
                    modelbinned_mat{s,ii} = modelmd_binned;
                else
                    modelmd_binned = modelbinned_mat{s,ii};
                end
                modelmd_binned_shuff = ARC_binAndTransform_shuffcoarse(modelmd_binned);
            else
                % Numerical control
                % 'Run this...'
                [modelmd_binned,modelmd_binned_mat] = ARC_binAndTransform_numctrl(modelmd2, behav_ratings_, binz, numctrl_bin);
                % modelmd_binned_shuff = ARC_binAndTransform_shuffcoarse_old(modelmd_binned);
                [~,modelmd_binned_shuff] = ARC_binAndTransform_shuffcoarse(modelmd2, behav_ratings_, binz, numctrl_bin,nshuff);
                modelmd_binned_shuff = permute(modelmd_binned_shuff,[3 1 2]);
            end

            if zscorer
                modelmd_binned = zscore(modelmd_binned,[],2);
            end
            modelmd_corrcoef = corrcoef(modelmd_binned);

            if sz_ctrl
                % Size control of ROIs
                [modelmd_corrcoef] = ARC_binAndTransform_sz(modelmd2, behav_ratings_, binz, [ms1 ms2],100);
                modelmd_binned_shuff = ARC_binAndTransform_shuffcoarse(modelmd_binned);
            end

            if intens_reg
                modelmd_binned_int = ARC_binAndTransform(modelmd2, behav_int(group_vec), binz, [min(behav_int) max(behav_int)]);
                % modelmd_binned_int  = zscore(modelmd_binned_int ,[],2);
                modelmd_corrcoef_int  = corrcoef(modelmd_binned_int);
                modelmd_binned_shuff_int = ARC_binAndTransform_shuffcoarse(modelmd_binned_int);
            end

            % Contruction of theoretical RSMs
            val_sc = linspace(-1,1,binz);
            sal_sc = abs(val_sc);
            valp = val_sc;
            valp(1:binzpart1)=0;
            valn = val_sc;
            valn(binzpart2:end)=0;


            val_mat = 1-abs(val_sc-val_sc');
            valp_mat = 1-abs(valp-valp');
            valn_mat = 1-abs(valn-valn');
            sal_mat = 1-abs(sal_sc-sal_sc');
            if run_sqdist
                valp_mat = 1-abs(valp-valp').^2;
                valn_mat = 1-abs(valn-valn').^2;
            end

            utl_mask1 = logical(blkdiag(zeros(binz-binzpart1),ones(binzpart1)));
            utl_mask2 = logical(triu(ones(length(val_mat)),1)); % All possible odors
            utl_mask = and(utl_mask1,utl_mask2);
            utl_mask_blk = flipud(utl_mask1);
            utl_mask_blk( binzpart1 , binzpart1 )=false;

            if noblock
                utl_mask_blk = utl_mask2;
            end

            if ~valsep % Basic valence and salience RSA
                if  intens_reg
                    % Intensity control
                    des_x = [val_mat(utl_mask2) sal_mat(utl_mask2) modelmd_corrcoef_int(utl_mask2)];
                elseif sniff_ctrl
                    fprintf('sniff ctrl in progress...')
                    des_x = [val_mat(utl_mask2) sal_mat(utl_mask2) sniff_corr{s}(utl_mask2)];
                    des_x_shuff = des_x;
                else
                    des_x = [val_mat(utl_mask2) sal_mat(utl_mask2) ];
                    des_x_shuff = des_x;
                end

                wt_dist = zeros(nshuff,size(des_x,2)+1);
                for zz = 1:nshuff % Bootstrap
                    modelmd_binned = squeeze(modelmd_binned_mat(:,:,zz));
                    if zscorer
                        modelmd_binned = zscore(modelmd_binned,[],2);
                    end
                    modelmd_corrcoef = corrcoef(modelmd_binned);
                    [wt_dist(zz,:) , ~] = ARC_multicomputeWeights_tsc(des_x,modelmd_corrcoef(utl_mask2));
                end
                wt = mean(wt_dist);

                wt_dist = zeros(nshuff,size(des_x,2)+1);
                for zz = 1:nshuff % Shuffle
                    modelmd_binneds = squeeze(modelmd_binned_shuff(zz,:,:));
                    if zscorer
                        modelmd_binneds = zscore(modelmd_binneds,[],2);
                    end
                    modelmd_corrcoef2 = corrcoef(modelmd_binneds);

                    if  intens_reg
                        int_shuff =  squeeze(modelmd_binned_shuff_int(zz,:,:));
                        if zscorer
                            int_shuff  = zscore( int_shuff ,[],2);
                        end
                        int_shuff_mat  = corrcoef( int_shuff );
                        des_x_shuff = [des_x(:,1:2) int_shuff_mat(utl_mask2)];
                    end
                    [wt_dist(zz,:), ~] = ARC_multicomputeWeights_tsc(des_x_shuff,modelmd_corrcoef2(utl_mask2));
                end
                rsa_Pps(s,ii,:,1:2) = wt_dist(:,2:3);
                rsa_Pp(s,ii,1) = 1-invprctile(wt_dist(:,2),wt(2))/100;
                rsa_Pp(s,ii,2) = 1-invprctile(wt_dist(:,3),wt(3))/100;
            else
                if intens_reg
                    des_x = [valp_mat(utl_mask_blk) valn_mat(utl_mask_blk)...
                        modelmd_corrcoef_int(utl_mask_blk) ];
                elseif sniff_ctrl
                    des_x = [valp_mat(utl_mask_blk) valn_mat(utl_mask_blk) sniff_corr{s}(utl_mask_blk)];
                    des_x_shuff = des_x;
                else

                    des_x = [valp_mat(utl_mask_blk) valn_mat(utl_mask_blk)];
                    des_x_shuff = des_x;
                end

                wt_dist = zeros(nshuff,size(des_x,2)+1);
                for zz = 1:nshuff % Bootstrap (num bootstrap = num shuffle)
                    modelmd_binned = squeeze(modelmd_binned_mat(:,:,zz));
                    if zscorer
                        modelmd_binned = zscore(modelmd_binned,[],2);
                    end
                    modelmd_corrcoef = corrcoef(modelmd_binned);
                    [wt_dist(zz,:) , ~] = ARC_multicomputeWeights_tsc( des_x,modelmd_corrcoef(utl_mask_blk));
                end
                wt = mean(wt_dist);
                % if (ii==3)
                %     'error'
                % end

                wt_dist = zeros(nshuff,size(des_x,2)+1);
                for zz = 1:nshuff % shuffle
                    modelmd_binneds = squeeze(modelmd_binned_shuff(zz,:,:));
                    if  zscorer
                        modelmd_binneds = zscore(modelmd_binneds,[],2);
                    end
                    modelmd_corrcoef2 = corrcoef(modelmd_binneds);

                    if  intens_reg
                        int_shuff =  squeeze(modelmd_binned_shuff_int(zz,:,:));
                        if zscorer
                            int_shuff  = zscore( int_shuff ,[],2);
                        end
                        int_shuff_mat  = corrcoef( int_shuff );
                        des_x_shuff = [des_x(:,1:2) int_shuff_mat(utl_mask_blk)];
                    end

                    [wt_dist(zz,:), ~] = ARC_multicomputeWeights_tsc( des_x_shuff,modelmd_corrcoef2(utl_mask_blk));
                end

                rsa_Pps(s,ii,:,1:2) = wt_dist(:,2:3);
                rsa_Pp(s,ii,1) = 1-invprctile(wt_dist(:,2),wt(2))/100;
                rsa_Pp(s,ii,2) = 1-invprctile(wt_dist(:,3),wt(3))/100;

            end
            rsa_P1wt(s,ii,:)  =wt(2:3);

            % Voxwise stuff
            % modelmd_binned = squeeze(mean(model_binned_mat,3));
            % t1 = zeros(nshuff,1);  t2 = zeros(nshuff,3); t3 = zeros(nshuff,1); t4 = zeros(size(modelmd_binned,1),2,nshuff);  t5 = zeros(size(modelmd_binned,1),2,nshuff);
            % for zz = 1:nshuff
            % modelmd_binned = squeeze(model_binned_mat(:,:,zz));
            % [ t1(zz),t4(:,:,zz), t2(zz,:),t5(:,:,zz),t3(zz)] = ARC_multicomputeWeights_tsc_voxwise(valp_mat, valn_mat, modelmd_binned);
            % end
            % rsa_Pcorr(s,ii) = mean(t1); t_scores = squeeze(mean(t4,3)); rsa_prop(s,ii,:) = ...
            %     mean(t2); w_scores = squeeze(mean(t5,3)); rsa_Pcorrp(s,ii) = mean(t3);
            % fprintf('Running snuffle')
            % modelmd_binned_mat = permute(modelmd_binned_shuff,[2 3 1]);
            %
            % if (ii==3)
            % 'error'
            % end
            if run_shuff
                fprintf('Running Shuff')
                modelmd_binned_mat = permute(modelmd_binned_shuff,[2 3 1]);
            end
            if run_voxwise
                [rsa_Pcorr(s,ii) , t_scores, rsa_prop(s,ii,:), w_scores, rsa_Pcorrp(s,ii), wt_mat] = ...
                    ARC_multicomputeWeights_tsc_voxwise_boot(valp_mat, valn_mat, modelmd_binned_mat,seq);
            else
                [rsa_Pcorr(s,ii) , t_scores, rsa_prop(s,ii,:), w_scores, rsa_Pcorrp(s,ii), wt_mat] = ...
                    ARC_multicomputeWeights_tsc_voxwise_boot(valp_mat, valn_mat, modelmd_binned_mat(:,:,1:2),seq);
            end

            w_score_mat{s,ii} = wt_mat ;
            wt_mat = permute(wt_mat,[2 1 3]);
            wt_mat_sort = reshape(wt_mat,[],2);
            % w_mat{ii} = cat(1,w_mat{ii},wt_mat_sort);
            w_mat{ii} = cat(1,w_mat{ii},w_scores);


            t_score_mat{s,ii} = t_scores;
            % w_score_mat{s,ii} = w_scores;

            df = length(t_scores);
            func = @(x) (1 - tcdf((x),df));
            p1_mat = arrayfun(func,t_scores(:,1));
            p2_mat = arrayfun(func,t_scores(:,2));

            thr_fdranat(s,1,ii) = tinv(1-fdr_benjhoc(p1_mat),df);
            thr_fdranat(s,2,ii) = tinv(1-fdr_benjhoc(p2_mat),df);

            col1Above = t_scores(:,1) >  thr_fdranat(s,1,ii)  & t_scores(:,2) <= thr_fdranat(s,2,ii); % Col 1 above thr, Col 2 not
            col2Above = t_scores(:,2) > thr_fdranat(s,2,ii) & t_scores(:,1) <=  thr_fdranat(s,1,ii) ; % Col 2 above thr, Col 1 not
            bothAbove = t_scores(:,1) >  thr_fdranat(s,1,ii)  & t_scores(:,2) > thr_fdranat(s,2,ii);  % Both cols above thr

            % Pos_population:
            if ~demomode
                rsa_pos_pop = nanmean(modelmd2( col1Above,:));
                rsa_neg_pop = nanmean(modelmd2( col2Above,:));

                % map_area(:,:,:,ii,1)  = unmasker(t_scores(:,1),logical(anatmasks(:,:,:,ii)));
                % map_area(:,:,:,ii,2) = unmasker(t_scores(:,2),logical(anatmasks(:,:,:,ii)));
            end


            %% Figures
            % % Make figure
            subplot(3,nanat,kk)
            hold on
            % ARC_scatterDensity(w_scores(:,1),w_scores(:,2))
            % xlabel('Val+ RSA beta')
            % ylabel('Val- RSA beta')
            % title(sprintf('S%02d: %s, r: %.2f',s,anat_names{ii},rsa_Pcorr(s,ii)))

            imagesc((val_sc),(val_sc),flipud(modelmd_corrcoef))
            colormap(ARC_customcmap(64))
            xticks(1:2*binz)
            yticks(1:2*binz)
            yticklabels(-1:0.5:0)
            xticklabels(0:0.5:1)
            kk = kk+1;
            colorbar
            xlabel('Pleasantness')
            ylabel('Pleasantness')
            title(anat_names{ii})
        end

        map_area(map_area==0)=nan;
        map_mat = squeeze(mean(map_area,4,"omitnan"));
        m1 = squeeze(map_mat(:,:,:,1));
        m2 = squeeze(map_mat(:,:,:,2));

        df1 = sum(~isnan(m1(:)));
        func = @(x) 2 * (1 - tcdf(abs(x),df1));
        p1_mat = arrayfun(func,m1(~isnan(m1)));

        df2 = sum(~isnan(m2(:)));
        func = @(x) 2 * (1 - tcdf(abs(x),df2));
        p2_mat = arrayfun(func,m2(~isnan(m2)));

        thr_fdr(s,1) = tinv(1-fdr_benjhoc(p1_mat),df1);
        thr_fdr(s,2) = tinv(1-fdr_benjhoc(p2_mat),df2);

        map_mat(isnan(map_mat))=0;
        mkdir(savepath)
        if ~demomode
            % write_reshaped_nifty(squeeze(map_mat(:,:,:,1)), savepath, false, fullfile(anatdir,maskfile), sprintf('ARC%02d_valp',s));
            % write_reshaped_nifty(squeeze(map_mat(:,:,:,2)), savepath, false, fullfile(anatdir,maskfile), sprintf('ARC%02d_valn',s));
        end
    end
end
savepath
savefig(fullfile(savepath,'imagescr'))
print(fullfile(savepath,'imagescr'),'-dpng')


%% Statistics and plotting
df = nchoosek(7,2);

% Average p-values
rsa_Pavg = zeros(nanat,2);
for ii=1:nanat
    for cond = 1:2
        meaneff = nanmean(rsa_P1wt(:,ii,cond));
        nulldist = nanmean(squeeze(rsa_Pps(:,ii,:,cond)),1);
        rsa_Pavg(ii,cond)= 1-invprctile(nulldist,meaneff)/100;
    end
end

% Average p-value voxelwise
rsa_Pavg_vox = zeros(nanat,1);
dfs = cellfun(@length,t_score_mat);
for ii=1:nanat
    meaneff = nanmean(rsa_Pcorr(:,ii));
    df = nanmean(dfs(:,ii));
    rsa_Pavg_vox(ii) = ARC_r2t( meaneff,df);
end

% Cross domain correlation
ARC_barplot_sig(rsa_P1wt,rsa_Pavg,true)
ylabel('RSA beta')
xticks(1:nanat)
xticklabels(anat_names)
if valsep
    legend({'Val+','Val-'})
else
    legend({'Valence','Salience'})
end
savefig(fullfile(savepath,'ARC_RSAwt'))
print(gcf,'-vector','-dsvg',[fullfile(savepath,'ARC_RSAwt'),'.svg']) % svg
print(fullfile(savepath,'ARC_RSAwt'),'-dpng')

% RSA correlation domains
figure()
hold on
bar(1:nanat,nanmean(rsa_Pcorr))
errorbar(1:nanat,nanmean(rsa_Pcorr),nanstd(rsa_Pcorr)/3,'.')
c_s = {'r','g','b'}; % Data dots for subjects
for jj = 1:3 % For bars for perceptual, chemical and combinations
    plot(1:nanat,rsa_Pcorr(jj,:),c_s{jj})
end
xticks(1:nanat)
xticklabels(anat_names);
ylabel('t score of correlation')
print(gcf,'-vector','-dsvg',[fullfile(savepath,'voxwise_similarity'),'.svg']) % svg
savefig(fullfile(savepath,'voxwise_similarity'))
print(fullfile(savepath,'voxwise_similarity'),'-dpng')

% RSA sig. voxels
ARC_barplot(rsa_prop)
xticks(1:nanat)
xticklabels(anat_names);
ylabel('Fraction voxels')
legend({'only Val+', 'only Val-', 'both'})
savefig(fullfile(savepath,'voxprop'))
print(fullfile(savepath,'voxprop'),'-dpng')
print(gcf,'-vector','-dsvg',[fullfile(savepath,'voxprop'),'.svg']) % svg


% Scatter density
figure('Position', [0.5 0.5 1420 240])
% clims = {[0 5000],[0 6000],[0 60000],[0 30000]};
clims = {[0 15],[0 20],[0 160],[0 80]};
hold on
kk = 0;
for ii = 1:nanat
    kk = kk+1;
    subplot(1,nanat,kk)
    hold on
    % ARC_scatterDensity(w_mat{ii}(:,1),w_mat{ii}(:,2))
    % histogram2(w_mat{ii}(:,1),w_mat{ii}(:,2),'DisplayStyle','tile','BinWidth',[0.05 0.05],'EdgeColor','none')

    w_mat = vertcat(w_score_mat{:,ii});

    % w_mat = reshape(w_mat,[],2);
    w_mat = squeeze(mean(w_mat,2));
    % vec = [-1:0.05:1];
    % vec(mode(discretize(w_mat, vec)))
    histogram2( w_mat(:,1), w_mat(:,2),'DisplayStyle','tile','BinWidth',[0.05 0.05],'EdgeColor','none')


    xlabel('Val+ RSA beta')
    ylabel('Val- RSA beta')
    xline(0)
    yline(0)
    % clim([0 6])
    xlim([-0.4 1])
    xticks([-0.4 0.3 1])
    ylim([-0.4 1])
    yticks([-0.4 0.3 1])
    colorbar
    % clim(clims{ii})
end
savefig(fullfile(savepath,'ARC_dens3'))
print(fullfile(savepath,'ARC_dens3'),'-dpng')
print(gcf,'-vector','-dsvg',[fullfile(savepath,'ARC_dens3'),'.svg']) % svg
save(fullfile(savepath,'ARC_wts'),'w_score_mat')
SFP_clearLargeVariables
save(fullfile(savepath,'ARC_RSA'))


%% Fig debug
% figure()
% hold on
% % ARC_scatterDensity(w_mat{ii}(:,1),w_mat{ii}(:,2))
% histogram2(ofc_means(:,1),ofc_means(:,2),'DisplayStyle','tile','BinWidth',[0.01 0.05],'EdgeColor','none')
% xlabel('Val+ RSA beta')
% ylabel('Val- RSA beta')
% xline(0)
% yline(0)
% % clim([0 6])
% xlim([-0.4 1])
% xticks([-0.4 0.3 1])
% ylim([-0.4 1])
% yticks([-0.4 0.3 1])


%% Fig percentiles
% load shuff_mean manually as w_score_mat2 and mainmat as w_score_mat
figure('Position', [0.5 0.5 1420 240])
hold on
kk = 0;
clims = {[0 20],[0 25],[0 200],[0 100]};
pvals= zeros(4,2);
for ii = 1:nanat
    kk = ii;
    subplot(3,nanat,kk)
    hold on
    % ARC_scatterDensity(w_mat{ii}(:,1),w_mat{ii}(:,2))
    cat1 = vertcat(w_score_mat{:,ii});
    cat2 = vertcat(w_score_mat2{:,ii});
    mat3 = ARC_inversePercentiles( cat1, cat2 ).*100;
    pvals(ii,1) = signrank( mat3(:,1), 0.1, 'tail', 'right');
    pvals(ii,2) = signrank( mat3(:,2), 0.1, 'tail', 'right');

    [~ ,pvals(ii,1)] =jbtest( mat3(:,1));
    [~,pvals(ii,2)] = jbtest( mat3(:,2));


    histogram2( mat3(:,1), mat3(:,2),'DisplayStyle','tile','BinWidth',[5 5],'EdgeColor','none')
    xlabel('Val+ RSA percentile')
    ylabel('Val- RSA percentile')
    colorbar
    xline(50)
    yline(50)
    clim(clims{ii})

    kk = ii+nanat;
    subplot(3,nanat,kk)
    histogram(mat3(:,1))
    xlabel('Val+ RSA percentile')
    ylabel('Freq.')

    kk = ii+2*nanat;
    subplot(3,nanat,kk)
    histogram(mat3(:,2))
    xlabel('Val- RSA percentile')
    ylabel('Freq.')

end

%% Percentile scores
% prctile_mat = cell(3,4);
% for ss = 1:3
% for ii = 1:nanat
%
%     % ARC_scatterDensity(w_mat{ii}(:,1),w_mat{ii}(:,2))
%     cat1 = (w_score_mat{ss,ii});
%     cat2 = (w_score_mat2{ss,ii});
%     prctile_mat{ss,ii} = ARC_inversePercentiles( cat1, cat2 ).*100;
%
% end
% end
pvals= zeros(4,2);
prctile_mat = cell(8,1);
for ii = 1:nanat
    % ARC_scatterDensity(w_mat{ii}(:,1),w_mat{ii}(:,2))
    cat1 = vertcat(w_score_mat{:,ii});
    cat2 = vertcat(w_score_mat2{:,ii});
    tmp = ARC_inversePercentiles( cat1, cat2 ).*100;
    prctile_mat{2*ii-1}= tmp(:,1);
    prctile_mat{2*ii}= tmp(:,2);

    % pvals(ii,1) = signrank( tmp(:,1), 50, 'tail', 'right');
    % pvals(ii,2) = signrank(tmp(:,2), 50, 'tail', 'right');

    % cat1_distp = cat1(:,:,1);
    % cat2_distp = cat2(:,:,1);
    % cat1_distn = cat1(:,:,2);
    % cat2_distn = cat2(:,:,2);
    % pvals(ii,1) = ranksum(cat1_distp(:),cat2_distp(:));
    % pvals(ii,2) = ranksum(cat1_distn(:),cat2_distn(:));

    % [~ ,pvals(ii,1)] =jbtest( tmp(:,1));
    % [~,pvals(ii,2)] = jbtest( tmp(:,2));

    [pvals(ii,1)] = signtest( (tmp(:,1)-50),0,'Tail','both');
    [pvals(ii,2)] = signtest( (tmp(:,2)-50),0,'Tail','both');
end


% %% Behavioral plots
% figure('Position',[0.5 0.5 960 250])
% v_ids = [2 2 2];
% fname = 'pleasantness';
% hold on
% bin_cent = [];
% for ss = 1:3
%     v_id = v_ids(ss);
%     if v_id<=18
%         subplot(1,3,ss)
%         behav_ratings = behav(ss).ratings(:,v_id);
%
%         behav_ratings = normalize(behav_ratings,'medianiqr');
%
%         [ms1,amin] = min(behav_ratings);
%         [ms2,amax] = max(behav_ratings);
%
%         edges  = linspace(ms1, ms2, 8);
%         cent = edges(1:end-1)+mean(diff(edges))/2;
%         histogram( behav_ratings,edges)
%         if ss==1
%             ylabel(fname)
%         end
%         title(sprintf('subject: %02d',ss))
%
%         xticks(round(cent,2))
%          val_sc = linspace(-1,1,7);
%
%          cellstr = cellfun(@(x) num2str(round(x,2)), num2cell(val_sc),'UniformOutput',false );
%         xticklabels( cellstr )
%         xtickangle(90)
%     end
% end
% savefig(fname)
% print(fname,'-dpng')
%
% % Check pdists
% thr1 = 1.98;
% thr2 = 1.98;
% figure('Position', [0.5 0.5 1280 240])
% edge = [0:0.5:20];
% for ss = 1:3
%     thr = tinv(0.90,size(t_score_mat1{s,3},1));
%     pospop1 = and(t_score_mat1{s,3}(:,1)>thr1,t_score_mat1{s,3}(:,2)<thr);
%     negpop1 = and(t_score_mat1{s,3}(:,2)>thr2,t_score_mat1{s,3}(:,1)<thr);
%
%     % thr1 = 1.65;
%     % thr2 = 1.65;
%     % pospop2 = and(t_score_mat2{s,3}(:,1)>thr1,t_score_mat2{s,3}(:,2)<thr);
%     % negpop2 = and(t_score_mat2{s,3}(:,2)>thr2,t_score_mat2{s,3}(:,1)<thr);
%
%     subplot(2,3,ss)
%     hold on
%     histogram(t_score_mat1{s,3}(pospop1,1),edge)
%     histogram(t_score_mat{s,3}(pospop1,1),edge)
%     title(sprintf('RSA val+ OFC Sub: %02d',ss))
%
%     subplot(2,3,ss+3)
%     hold on
%     histogram(t_score_mat1{s,3}(negpop1,2),edge)
%     histogram(t_score_mat{s,3}(negpop1,2),edge)
%     title(sprintf('RSA val- OFC Sub: %02d',ss))
% end
% savefig(fullfile(savepath,'distcompar'))
