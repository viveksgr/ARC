%% General Settings

% Main settings adjusted: val_sep, fmasker, control conditions
% shuffler

root = 'C:\Work\ARC\ARC';
maskfile =  'ARC3_anatgw.nii';
fmaskfile = 'ARC3_fanatgw3_pos.nii';
fmasker = false;
binz =7 ;
betas = -2:0.1:2;
if mod(binz,2)==0; binzpart1 = binz/2; binzpart2 = binzpart1+1; else; binzpart1 = (binz+1)/2 ; binzpart2 = binzpart1; end

anat_names = {'PC','AMY','OFC','VMPFC'};
anat_masks = {'rwPC.nii','rwAmygdala.nii','rwofc.nii','rwvmpfc.nii'};
nanat = length(anat_names);

num_cntrl = false;
sz_ctrl = false;
intens_reg = true;

valsep = true;
single_n = false; % Noisepool
single_c = true; % Cutoff from sign voxels
zscorer = true;
rangenormer =false;
noblock = false;

nodor = 160;
nshuff = 1000;
shuffler = false;
dirs = {fullfile(root,'\ARC01\mediation');
    fullfile(root,'\ARC02\mediation');
    fullfile(root,'\ARC03\mediation')};
behav = load(fullfile(root,'ARC','NEMO_perceptual2.mat'));
modelname = '\unmasked\control\temp';
savepath = fullfile(root,'RSA',modelname);
v_id = 2; % Index of vector for median splitting odors

% load(fullfile(statpath,'fir_cv.mat'))
fprintf('\n')

% Result matrics
rsa_P1 = zeros(3,nanat,2); % For t-scores in anatomical areas
rsa_P1wt = zeros(3,nanat,length(betas)); % For beta-wts in anatomical areas
rsa_Pp = zeros(3,nanat,2); %
rsa_Pps = zeros(3,nanat,nshuff,length(betas)-1); % Permutation test
rsa_Pcorr = zeros(3,nanat); % Correlation of beta weights
rsa_prop = zeros(3,nanat,3);  % Proportion of significant voxels

t_score_mat = cell(3,4);
w_mat = cell(nanat,1);
% Stats
thr_fdr = zeros(3,2);
thr_fdranat = zeros(3,2,nanat);

% Figure across subjects
figure1 = figure('OuterPosition',[297 183 1209 737]);
hold on
kk = 1;

for s = [1 2 3] % Subject
    fprintf('Subject: %02d\n',s)
    anatdir = fullfile(root,sprintf('ARC%02d',s),'single');
    % Construction of median split
    behav_ratings = behav.behav(s).ratings(:,v_id);

    % behav_ratings = behav_ratings(randperm(length(behav_ratings)));
    behav_int = behav.behav(s).ratings(:,1);
    behav_ratings = normalize(behav_ratings,'medianiqr');
    ms1 = min(behav_ratings);
    ms2 = max(behav_ratings);
    behav_int = normalize( behav_int,'medianiqr');
    md = 0;

    statpath = dirs{s};
    % Gray Matter, Functional and Anatomical Masks
    mask = (spm_read_vols(spm_vol(fullfile(statpath, maskfile)))); % Mask used to construct odor files
    mask(isnan(mask))=0;
    mask = logical(mask);
    fmask = (spm_read_vols(spm_vol(fullfile(statpath, fmaskfile)))); % Mask used to examine voxels in RSA
    fmask(isnan(fmask))=0;
    if ~fmasker
        fmask = fmask +0.1;
    end
    fmask = logical(fmask); % Only choose voxels with significant odor evoked activity
    fmask_1d = fmask(mask);
    marea = and(fmask,mask);

    if s <3
        anatpath = fullfile(root,sprintf('NEMO_%02d',s),'\imaging\nii\masks');
    else
        anatpath = fullfile(root,sprintf('NEMO_%02d',s+1),'\imaging\nii\masks');
    end

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

    if s==3; s2 = 4; else; s2 = s; end
    onsets = load(fullfile(anatdir,sprintf('conditions_NEMO%02d.mat',s2)),'onsets');
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

    if shuffler
        behav_ratings_ = behav_ratings_(randperm(length(behav_ratings_)));
    end

    for ii = 1:length(anat_names)

        fprintf('area:%02d\n',ii)
        modelmd_ = load(fullfile(anatdir,anat_names{ii},'TYPED_FITHRF_GLMDENOISE_RR.mat'),'modelmd','noisepool');
        modelmd = squeeze(modelmd_.modelmd);
        noisepool = modelmd_.noisepool;
        if single_c
            modelmd = modelmd(masks_set_cell{ii},:);
            noisepool = noisepool(masks_set_cell{ii});
        end

        if single_n
            modelmd2 = modelmd(~noisepool,:);
        else
            modelmd2 = modelmd;
        end

        fprintf('size:%02d\n',size(modelmd,1))

        [r1,~] = find(isnan(modelmd2));
        modelmd2(r1,:) = [];

        % if s==3
        %     'beep'
        % end

        % [M1_new, M2_new] = ARC_transformMatrix(behav_ratings_,  modelmd2, binz,[ms1 ms2]);
        if ~rangenormer
            if ~num_cntrl
                % Basic model
                modelmd_binned = ARC_binAndTransform(modelmd2, behav_ratings_, binz, [ms1 ms2]);
                % modelmd_binned_shuff = ARC_binAndTransform_shuff(modelmd2, behav_ratings_, binz, [ms1 ms2],nshuff);
                modelmd_binned_shuff = ARC_binAndTransform_shuffcoarse(modelmd_binned);
            else
                modelmd_binned = ARC_binAndTransform_numctrl(modelmd2, behav_ratings_, binz, [ms1 ms2]);
                modelmd_binned_shuff = ARC_binAndTransform_shuffcoarse(modelmd_binned);
            end
        else
            modelmd_binned = ARC_binAndTransformQuantiles(modelmd2, behav_ratings_, binz);
            % modelmd_binned_shuff = ARC_binAndTransformQuantiles_shuff(modelmd2, behav_ratings_, binz,nshuff);
            modelmd_binned_shuff = ARC_binAndTransform_shuffcoarse(modelmd_binned);
            edges_ = quantile(behav_ratings_, linspace(0, 1, binz+1));
            centers_ = edges_(1:end-1)+(edges_(2)-edges_(1))/2;
            centers_rn = 2*((centers_-ms1)./(ms2-ms1))-1;
            centers_rd = round(centers_rn,2);
        end
        if zscorer
            modelmd_binned = zscore(modelmd_binned,[],2);
        end
        modelmd_corrcoef = corrcoef(modelmd_binned);

        if sz_ctrl
            [modelmd_corrcoef] = ARC_binAndTransform_sz(modelmd2, behav_ratings_, binz, [ms1 ms2],100);
            modelmd_binned_shuff = ARC_binAndTransform_shuffcoarse(modelmd_binned);
        end

        if intens_reg
            modelmd_binned_int = ARC_binAndTransform(modelmd2, behav_int(group_vec), binz, [min(behav_int) max(behav_int)]);
            % modelmd_binned_int  = zscore(modelmd_binned_int ,[],2);
            modelmd_corrcoef_int  = corrcoef(modelmd_binned_int);
            modelmd_binned_shuff_int = ARC_binAndTransform_shuffcoarse(modelmd_binned_int);
        end

        des_x = zeros(nchoosek(binz,2),length(betas));
        utl_mask1 = logical(blkdiag(zeros(binz-binzpart1),ones(binzpart1)));
        utl_mask2 = logical(triu(ones(binz),1)); % All possible odors

        val_sc = linspace(-1,1,binz);
        betas = -2:0.1:2;
        des_x = zeros(nchoosek(binz,2),length(betas));
        val_par = zeros(binz,length(betas));
        figure()
        hold on
        for tt = 1:length(betas)-1
            subplot(5,8,tt)
            hold on
            scale = betas(tt)*ones(binz,1);
            scale(ceil((binz-1)/2)+1:end)=1;
            val_par(:,tt)=scale.*val_sc';
            imagemat = 1-abs(val_par(:,tt)-val_par(:,tt)');
            des_x(:,tt)= imagemat(utl_mask2);
            imagesc(linspace(-1,1,binz),linspace(-1,1,binz),imagemat )
            title(sprintf('slope: %.2f',betas(tt)))
            axis tight
        end
        % imagesc(val_par)

        [wt, t_scores] = ARC_multicomputeWeights_onesc( des_x,modelmd_corrcoef(utl_mask2));
        
        % y = regressmeout(modelmd_corrcoef(utl_mask2)',modelmd_corrcoef_int(utl_mask2)')';
        % [wt, t_scores] = ARC_multicomputeWeights_tsc(  [val_mat(utl_mask2) sal_mat(utl_mask2) ],y);

        wt_dist = zeros(nshuff,size(des_x,2)+1);
        for zz = 1:nshuff
            modelmd_binneds = squeeze(modelmd_binned_shuff(zz,:,:));
            if zscorer
                modelmd_binneds = zscore(modelmd_binneds,[],2);
            end
            modelmd_corrcoef2 = corrcoef(modelmd_binneds);
            [wt_dist(zz,:), ~] = ARC_multicomputeWeights_onesc( des_x,modelmd_corrcoef2(utl_mask2));
        end

        rsa_Pps(s,ii,:,1:end) = wt_dist(:,2:end-1);
        % rsa_Pp(s,ii,1) = 1-invprctile(wt_dist(:,2),wt(2))/100;
        % rsa_Pp(s,ii,2) = 1-invprctile(wt_dist(:,3),wt(3))/100;

        rsa_P1wt(s,ii,:)  = wt(2:end);
    end
end
ARC_plotMatrixLines(rsa_P1wt,anat_names)

% figure()
% plot(squeeze(mean(rsa_P1wt))')

% ARC_barplot(rsa_P1wt)
% df = nchoosek(7,2);
% p_values_3d = ARC_RSA_pvals(rsa_P1, rsa_P1wt, df)
% rsa_Pp
%
% % pvalue avg
rsa_Pavg = zeros(nanat,length(betas)-1);
for ii=1:nanat
    for cond = 1:length(betas)-1
        meaneff = mean(rsa_P1wt(:,ii,cond));
        nulldist = mean(squeeze(rsa_Pps(:,ii,:,cond)),1);
        rsa_Pavg(ii,cond)= 1-invprctile(nulldist,meaneff)/100;
    end
end

figure()
mean_eff = squeeze(mean(rsa_P1wt))';
plot(betas(1:end-1),mean_eff(1:end-1,:))
hold on

sig_t = rsa_Pavg>0.05;
mean_eff(sig_t)=nan;
plot(betas(1:end-1),mean_eff(1:end-1,:),'Linestyle',"-",'LineWidth',2)

% % rsa_Pavg
% rsa_Pavg_vox = zeros(nanat,1);
% dfs = cellfun(@length,t_score_mat);
% for ii=1:nanat
%         meaneff = mean(rsa_Pcorr(:,ii))
%         df = mean(dfs(:,ii));
%        rsa_Pavg_vox(ii) = ARC_r2t( meaneff,df);
%
% end
% % Cross domain correlation
%
% mkdir(savepath)
% % RSA beta plot
% ARC_barplot_sig(rsa_P1wt,rsa_Pavg)
% % ARC_barplot_sig(rsa_P1wt,p_values_3d)
% ylabel('RSA beta')
% xticks(1:nanat)
% xticklabels(anat_names)
% if valsep
%     legend({'Val+','Val-'})
% else
%     legend({'Valence','Salience'})
% end
% savefig(fullfile(savepath,'ARC_RSAwt'))
% print(gcf,'-vector','-dsvg',[fullfile(savepath,'ARC_RSAwt'),'.svg']) % svg
% print(fullfile(savepath,'ARC_RSAwt'),'-dpng')
%
% % RSA correlation domains
% figure()
% hold on
% bar(1:nanat,mean(rsa_Pcorr))
% errorbar(1:nanat,mean(rsa_Pcorr),std(rsa_Pcorr)/3,'.')
% c_s = {'r','g','b'}; % Data dots for subjects
% for jj = 1:3 % For bars for perceptual, chemical and combinations
%     plot(1:nanat,rsa_Pcorr(jj,:),c_s{jj})
% end
% xticks(1:nanat)
% xticklabels(anat_names);
% ylabel('t score of correlation')
% print(gcf,'-vector','-dsvg',[fullfile(savepath,'voxwise_similarity'),'.svg']) % svg
% savefig(fullfile(savepath,'voxwise_similarity'))
% print(fullfile(savepath,'voxwise_similarity'),'-dpng')
%
% % RSA sig. voxels
% figure()
% ARC_barplot(rsa_prop)
% xticks(1:nanat)
% xticklabels(anat_names);
% ylabel('Fraction voxels')
% legend({'only Val+', 'only Val-', 'both'})
% savefig(fullfile(savepath,'voxprop'))
% print(fullfile(savepath,'voxprop'),'-dpng')
% print(gcf,'-vector','-dsvg',[fullfile(savepath,'voxprop'),'.svg']) % svg
%
% % Scatter density
% figure()
% hold on
% kk = 0;
% for ii = 1:4
%     kk = kk+1;
%     subplot(1,4,kk)
%     hold on
%     ARC_scatterDensity(w_mat{ii}(:,1),w_mat{ii}(:,2))
%     xlabel('Val+ RSA beta')
%     ylabel('Val- RSA beta')
%     clim([0 1.5])
%     xlim([-0.4 1])
%     xticks([-0.4 0.3 1])
%     ylim([-0.4 1])
%     yticks([-0.4 0.3 1])
% end
% savefig(fullfile(savepath,'ARC_dens'))
% print(fullfile(savepath,'ARC_dens'),'-dpng')
% print(gcf,'-vector','-dsvg',[fullfile(savepath,'ARC_dens'),'.svg']) % svg
% SFP_clearLargeVariables
% save(fullfile(savepath,'ARC_RSA'))
%
%
% % %% Visualize anat RDS
% % imger = false;
% % if imger
% %     [bsort,argsort ] = sort(behav_ratings_);
% %     ii = 2;
% %     % -- Find M_anat_amy
% %     M_anat_amy = M_anat(argsort,argsort);
% %     figure()
% %     subplot(1,2,1);
% %     imagesc(M_anat(argsort,argsort));
% %     ii = 3;
% %     % -- Repeat for ofc
% %     M_anat_ofc = M_anat(argsort,argsort);
% %
% %     % Smooth 2D
% %     [~,histedge] = histcounts(bsort,40);
% %     histparam = find_nearest(bsort,histedge);
% %
% %     M_tar = M_anat_amy;
% %     M_sm_amy= [];
% %     for ii = 1:length(histparam)-1
% %         for jj = 1:length(histparam)-1
% %             temp = M_tar(histparam(ii):histparam(ii+1),histparam(jj):histparam(jj+1));
% %             if ~isempty(temp)
% %                 M_sm_amy(ii,jj) = mean(temp(:));
% %             end
% %         end
% %     end
% %
% %     % plot
% %     figure()
% %     subplot(1,2,1);
% %     imagesc(M_sm_amy);
% %     caxis([0.97 1.03])
% %     colorbar
% %     title('Amy')
% %     subplot(1,2,2);
% %     imagesc(M_sm_ofc);
% %     caxis([0.97 1.03])
% %     colorbar
% %     title('ofc')
% %
% %     % Valence RDM
% %     idx2 = std(M_val_img,1)==0;
% %     M_red = M_val_img(~idx2,~idx2);
% %     [idx] = kmeans(M_red,5,'Distance','correlation');
% %     [kidx] = kmeans(behav_ratings_,4);
% %     [~,argsort] =sort(kidx);
% %     imagesc( M_val(argsort,argsort))
% %     colormap(swampsunset)
% %     idx_exp = 1:1:length(M_val_img);
% %     idx_exp_km = idx_exp(idx);
% %     idx_exp_nkm = idx_exp(idx2);
% %     idx_idx = [idx_exp_km, idx_exp_nkm];
% %     M_val_img = M_val_img(idx_idx, idx_idx);
% %     imagesc(M_val_img)
% % end
% %
% % %% RDMs
% % figure()
% % hold on
% % subplot(2,4,1)
% % plot(val_sc,val_sc)
% % ylabel('Valence')
% %
% % subplot(2,4,2)
% % plot(val_sc,sal_sc)
% % ylabel('Salience')
% %
% % subplot(2,4,3)
% % plot(val_sc,valp)
% % ylabel('Val+')
% %
% % subplot(2,4,4)
% % plot(val_sc,valn)
% % ylabel('Val-')
% %
% % subplot(2,4,5)
% % val_mat = 1-abs(val_sc-val_sc');
% % valp_mat = 1-abs(valp-valp');
% % valn_mat = 1-abs(valn-valn');
% % sal_mat = 1-abs(sal_sc-sal_sc');
% % imagesc(val_sc,val_sc,flipud(val_mat))
% % colormap(ARC_customcmap(64))
% % yticks([-1 -0.5 0 0.5 1])
% % yticklabels(-[-1 -0.5 0 0.5 1])
% % xlabel('Pleasantness')
% % ylabel('Pleasantness')
% %
% % subplot(2,4,6)
% % imagesc(val_sc,val_sc,flipud(sal_mat))
% % colormap(ARC_customcmap(64))
% % yticks([-1 -0.5 0 0.5 1])
% % yticklabels(-[-1 -0.5 0 0.5 1])
% % xlabel('Pleasantness')
% %
% % subplot(2,4,7)
% % imagesc(val_sc,val_sc,flipud(valp_mat))
% % colormap(ARC_customcmap(64))
% % yticks([-1 -0.5 0 0.5 1])
% % yticklabels(-[-1 -0.5 0 0.5 1])
% % xlabel('Pleasantness')
% %
% % subplot(2,4,8)
% % imagesc(val_sc,val_sc,flipud(valn_mat))
% % colormap(ARC_customcmap(64))
% % yticks([-1 -0.5 0 0.5 1])
% % yticklabels(-[-1 -0.5 0 0.5 1])
% % xlabel('Pleasantness')
% % colorbar
% %
% % %% Val Sal
% % figure()
% % for s = [1 2 3] % Subject
% %     fprintf('Subject: %02d\n',s)
% %     anatdir = fullfile(root,sprintf('ARC%02d',s),'single');
% %     % Construction of median split
% %     behav_ratings = behav.behav(s).ratings(:,v_id);
% %     behav_int = behav.behav(s).ratings(:,1);
% %     % behav_ratings = vs_normalizer(behav_ratings);
% %
% %     % Valence and Salience distributions
% %     hold on
% %     subplot(2,3,s)
% %     histogram(behav_ratings,20)
% %
% %     subplot(2,3,s+3)
% %     histogram(abs(behav_ratings),20)
% % end
% %
%
% % Avg behavior
% cid = [behav.behav(1).cid; behav.behav(2).cid; behav.behav(4).cid];
% bv = [behav.behav(1).ratings(:,2); behav.behav(2).ratings(:,2); behav.behav(4).ratings(:,2)];
% cid_g = findgroups(cid);
% val_avg = splitapply(@nanmean,bv,cid_g);
