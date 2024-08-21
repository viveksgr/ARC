%% General Settings
root = 'C:\Work\ARC\ARC';
maskfile =  'ARC3_anatgw.nii';
fmaskfile = 'ARC3_fanatgw3_pos.nii';
fmasker = true;
binz =7 ;
if mod(binz,2)==0; binzpart1 = binz/2; binzpart2 = binzpart1+1; else; binzpart1 = (binz+1)/2 ; binzpart2 = binzpart1; end

beta_wter = true;
anat_names = {'PC','AMY','OFC','VMPFC'};
anat_masks = {'rwPC.nii','rwAmygdala.nii','rwofc.nii','rwvmpfc.nii'};
% anat_names = {'PC','AMY','Thal','frontal_sal','frontal_val'};
% anat_masks = {'rwPC.nii','rwAmygdala.nii','rwThal.nii','frontal_sal.nii','frontal_val.nii'};
nanat = length(anat_names);
medianize_behav = true;

rangenormer = false;
num_cntrl = false;
sz_ctrl = false;
raw_RSA = false;
valsep =true;
intens_reg = true;
    
% sess_l = cat(3,nchoosek([1 2 3],2),nchoosek([2 3 4],2),nchoosek([2 3
% 4],2),nchoosek([2 3 4],2)); % For sesswise
% load('C:\Data\NEMO\swampsunset.mat');

single_trial = true;
single_z = true; % Zscore
coarser_ = false; % Only binary estimation of valence RSM
single_n = false; % Noisepool
single_c = true; % Cutoff from sign voxels
zscorer = false;
increm = false; % Compare rating to previous
partial_corr = false;

nodor = 160;


sess_l = repmat([0],1,2,3);
dirs = {fullfile(root,'\ARC01\mediation');
    fullfile(root,'\ARC02\mediation');
    fullfile(root,'\ARC03\mediation')};
behav = load(fullfile(root,'ARC','NEMO_perceptual2.mat'));
% modelname = fullfile('mainval',sprintf('bin%01d',binz));
modelname = 'controls\int_domainp';
savepath = fullfile(root,'RSA',modelname);
v_id = 2; % Index of vector for median splitting odors

% load(fullfile(statpath,'fir_cv.mat'))
fprintf('\n')

rsa_P1 = zeros(3,nanat,2);
rsa_P1wt = zeros(3,nanat,2);
rsa_Pcorr = zeros(3,nanat);
rsa_prop = zeros(3,nanat,3);

rs = zeros(3,1);
hold on
% Subject - index
t_score_mat = cell(3,4);
kk = 1;
figure1 = figure('OuterPosition',[297 183 1209 737]);
hold on
w_mat = cell(4,1);
thr_fdr = zeros(3,2);
thr_fdranat = zeros(3,2,nanat);

for s = [1 2 3] % Subject
    fprintf('Subject: %02d\n',s)
    anatdir = fullfile(root,sprintf('ARC%02d',s),'single');
    % Construction of median split
    behav_ratings = behav.behav(s).ratings(:,v_id);
    behav_int = behav.behav(s).ratings(:,1);
    % behav_ratings = vs_normalizer(behav_ratings);
    if medianize_behav
        behav_ratings = normalize(behav_ratings,'medianiqr');
        ms1 = min(behav_ratings);
        ms2 = max(behav_ratings);
         behav_int = normalize( behav_int,'medianiqr');
    else
        ms1 = -1;
        ms2 = 1;
    end
    % behav_ratings = zscore(behav_ratings);
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

    masks_set(isnan(masks_set))=0;
    linux_config = false;
    warning('off','all')
    map_area = zeros([size(anat_cell{s}) 2]);
    %% Representational connectivity
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

        if zscorer
            modelmd2 = zscore(modelmd2,[],2);
        end
        % [M1_new, M2_new] = ARC_transformMatrix(behav_ratings_,  modelmd2, binz,[ms1 ms2]);


        if ~rangenormer
            if ~num_cntrl
                modelmd_binned = ARC_binAndTransform(modelmd2, behav_ratings_, binz, [ms1 ms2]);

            else
                modelmd_binned = ARC_binAndTransform_sz(modelmd2, behav_ratings_, binz, [ms1 ms2]);
            end
        else
            modelmd_binned = ARC_binAndTransformQuantiles(modelmd2, behav_ratings_, binz);
            edges_ = quantile(behav_ratings_, linspace(0, 1, binz));
            centers_ = edges_(1:end-1)+(edges_(2)-edges_(1))/2;
            centers_rn = 2*((centers_-ms1)./(ms2-ms1))-1;
            centers_rd = round(centers_rn,2);
        end
        % modelmd_binned = zscore(modelmd_binned,[],2);
        modelmd_corrcoef = corrcoef(modelmd_binned);

        if sz_ctrl
            [modelmd_corrcoef] = ARC_binAndTransform_sz(modelmd2, behav_ratings_, binz, [ms1 ms2],100);
        end

        if intens_reg
            modelmd_binned_int = ARC_binAndTransform(modelmd2, behav_int(group_vec), binz, [min(behav_int) max(behav_int)]);
            % modelmd_binned_int  = zscore(modelmd_binned_int ,[],2);
            modelmd_corrcoef_int  = corrcoef(modelmd_binned_int);
        end

        % if and(s==2,ii==2)
        %     'dil ki surkh deewaron pe...'
        % end

        val_sc = linspace(-1,1,binz);
        if rangenormer;  val_sc = centers_rn; end

        sal_sc = abs(val_sc);
        valp = val_sc;
        valp(1:binzpart1)=0;
        valn = val_sc;
        valn(binzpart2:end)=0;

        val_mat = 1-abs(val_sc-val_sc');
        valp_mat = 1-abs(valp-valp');
        valn_mat = 1-abs(valn-valn');
        sal_mat = 1-abs(sal_sc-sal_sc');

        % % Non-linear
        % % val_sc_sq = sign(val_sc).*val_sc.^2;
        % val_sc_sq = sign(val_sc).*val_sc.^2;
        % valsq_mat = 1-abs(val_sc_sq-val_sc_sq');
        % valp = val_sc_sq;
        % valp(1:binz)=0;
        % valn = val_sc_sq;
        % valn(binz+1:end)=0;

        % valp = exp(val_sc);
        % valn = exp(-val_sc);
        % valp_mat = 1-abs(valp-valp');
        % valn_mat = 1-abs(valn-valn');


        utl_mask1 = logical(blkdiag(zeros(binz-binzpart1),ones(binzpart1)));
        % utl_mask1 = logical(blkdiag(zeros(binz-(binzpart2-1)),ones(binzpart2-1)));
        utl_mask2 = logical(triu(ones(length(val_mat)),1)); % All possible odors
        utl_mask = and(utl_mask1,utl_mask2);
        utl_mask_blk = flipud(utl_mask1);

        if ~valsep
            [wt, t_scores] = ARC_multicomputeWeights_tsc( [val_mat(utl_mask2) sal_mat(utl_mask2) ],modelmd_corrcoef(utl_mask2));
            if  intens_reg
                [wt, t_scores] = ARC_multicomputeWeights_tsc( [val_mat(utl_mask2) sal_mat(utl_mask2) modelmd_corrcoef_int(utl_mask2)],modelmd_corrcoef(utl_mask2));
            end
        else
            [wt, t_scores] = ARC_multicomputeWeights_tsc( [valp_mat(utl_mask_blk) valn_mat(utl_mask_blk)],modelmd_corrcoef(utl_mask_blk));
            if  intens_reg
                [wt, t_scores] = ARC_multicomputeWeights_tsc( [valp_mat(utl_mask_blk) valn_mat(utl_mask_blk)...
                    modelmd_corrcoef_int(utl_mask_blk) ],modelmd_corrcoef(utl_mask_blk));
            end
        end
        rsa_P1(s,ii,:)  = t_scores(2:3);
        rsa_P1wt(s,ii,:)  =wt(2:3);


        % utl_mask = logical(triu(ones(length(val_mat)),1)); % All possible odors
        % rsa_P1(s,ii,1) = fastcorr(modelmd_corrcoef,val_mat);

        [rsa_Pcorr(s,ii),t_scores, rsa_prop(s,ii,:),w_scores] = ARC_multicomputeWeights_tsc_voxwise(valp_mat, valn_mat, modelmd_binned);
        w_mat{ii} = cat(1,w_mat{ii},w_scores);
        t_score_mat{s,ii} = t_scores;
        % Calculate conditions
        % thr = tinv(0.975,size(t_scores,1));

        df = length(t_scores);
        func = @(x) 2 * (1 - tcdf(abs(x),df));   
        p1_mat = arrayfun(func,t_scores(:,1));
        p2_mat = arrayfun(func,t_scores(:,2));

        thr_fdranat(s,1,ii) = tinv(1-fdr_benjhoc(p1_mat),df);
        thr_fdranat(s,2,ii) = tinv(1-fdr_benjhoc(p2_mat),df);
        
        col1Above = t_scores(:,1) >  thr_fdranat(s,1,ii)  & t_scores(:,2) <= thr_fdranat(s,2,ii); % Col 1 above thr, Col 2 not
        col2Above = t_scores(:,2) > thr_fdranat(s,2,ii) & t_scores(:,1) <=  thr_fdranat(s,1,ii) ; % Col 2 above thr, Col 1 not
        bothAbove = t_scores(:,1) >  thr_fdranat(s,1,ii)  & t_scores(:,2) > thr_fdranat(s,2,ii);  % Both cols above thr

        % Pos_population:
        rsa_pos_pop = nanmean(modelmd2( col1Above,:));
        rsa_neg_pop= nanmean(modelmd2( col2Above,:));

        % rsa_Pcorr(s,ii) = corr(rsa_pos_pop', rsa_neg_pop');
        % Make figure
        subplot(3,4,kk)
        hold on

        ARC_scatterDensity(w_scores(:,1),w_scores(:,2))
        xlabel('Val+ RSA beta')
        ylabel('Val- RSA beta')
        title(sprintf('S%02d: %s, r: %.2f',s,anat_names{ii},rsa_Pcorr(s,ii)))

        % % imagesc((val_sc),(val_sc),flipud(modelmd_corrcoef-mean(modelmd_corrcoef,'all')))
        % figure()
        % modelmd_corrcoef(logical(eye(size(modelmd_corrcoef))))=nan;
        % imagesc((val_sc),(val_sc),(modelmd_corrcoef))
        % colormap(ARC_customcmap(64))
        % yticks([-1,0,1])
        % yticklabels({1,0,-1})
        % % imagesc(M_mat)
        %
        % %
        if ~rangenormer
            xticks(1:2*binz)
            yticks(1:2*binz)
            yticklabels(-1:0.5:0)
            xticklabels(0:0.5:1)
        end
        %
        kk = kk+1;
        % colorbar
        % axis tight
        % clim([-inf  0.7])
        % xlabel('Pleasantness')
        % ylabel('Pleasantness')
        % title(anat_names{ii})

        map_area(:,:,:,ii,1)  = unmasker(t_scores(:,1),logical(anatmasks(:,:,:,ii)));
        map_area(:,:,:,ii,2) = unmasker(t_scores(:,2),logical(anatmasks(:,:,:,ii)));
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
    write_reshaped_nifty(squeeze(map_mat(:,:,:,1)), savepath, false, fullfile(anatdir,maskfile), sprintf('ARC%02d_valp',s));
    write_reshaped_nifty(squeeze(map_mat(:,:,:,2)), savepath, false, fullfile(anatdir,maskfile), sprintf('ARC%02d_valn',s));
  
end
savefig(fullfile(savepath,'imagescr'))
print(fullfile(savepath,'imagescr'),'-dpng')
df = nchoosek(7,2);
p_values_3d = ARC_RSA_pvals(rsa_P1, rsa_P1wt, df)

mkdir(savepath)



ARC_barplot(rsa_P1)
ylabel('RSA t')
if valsep
    legend({'Val+','Val-'})
else
    legend({'Valence','Salience'})
end
clear modelmd_ modelmd modelmd2 S1_omat_vals S2_omat_vals unity M_anat M_sal M_val
savefig(fullfile(savepath,'ARC_RSA'))
save(fullfile(savepath,'ARC_RSA'))
print(fullfile(savepath,'ARC_RSA'),'-dpng')


ARC_barplot(rsa_P1wt)
ylabel('RSA beta')
if valsep
    legend({'Val+','Val-'})
else
    legend({'Valence','Salience'})
end
savefig(fullfile(savepath,'ARC_RSAwt'))
print(fullfile(savepath,'ARC_RSAwt'),'-dpng')


figure()
hold on
bar(1:nanat,mean(rsa_Pcorr))
errorbar(1:nanat,mean(rsa_Pcorr),std(rsa_Pcorr)/3,'.')
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

figure()
ARC_barplot(rsa_prop)
xticks(1:nanat)
xticklabels(anat_names);
ylabel('Fraction voxels')
legend({'only Val+', 'only Val-', 'both'})
savefig(fullfile(savepath,'voxprop'))
print(fullfile(savepath,'voxprop'),'-dpng')
print(gcf,'-vector','-dsvg',[fullfile(savepath,'voxprop'),'.svg']) % svg

% figure()
% hold on
% for ii = 1:nanat
%     subplot(2,nanat,ii)
%     plot(squeeze(rsa_pos_pop(:,ii,:))')
%     subplot(2,nanat,ii+nanat)
%       plot(squeeze(rsa_neg_pop(:,ii,:))')
% end

figure()
hold on
kk = 0;
for ii = 1:4
    kk = kk+1;
    subplot(1,4,kk)
    hold on
    ARC_scatterDensity(w_mat{ii}(:,1),w_mat{ii}(:,2))
    xlabel('Val+ RSA beta')
    ylabel('Val- RSA beta')
    clim([0 1.5])
    xlim([-0.4 1])
    xticks([-0.4 0.3 1])
    ylim([-0.4 1])
    yticks([-0.4 0.3 1])
end
savefig(fullfile(savepath,'ARC_dens'))
print(fullfile(savepath,'ARC_dens'),'-dpng')
print(gcf,'-vector','-dsvg',[fullfile(savepath,'ARC_dens'),'.svg']) % svg

%% Visualize anat RDS
imger = false;
if imger
    [bsort,argsort ] = sort(behav_ratings_);
    ii = 2;
    % -- Find M_anat_amy
    M_anat_amy = M_anat(argsort,argsort);
    figure()
    subplot(1,2,1);
    imagesc(M_anat(argsort,argsort));
    ii = 3;
    % -- Repeat for ofc
    M_anat_ofc = M_anat(argsort,argsort);

    % Smooth 2D
    [~,histedge] = histcounts(bsort,40);
    histparam = find_nearest(bsort,histedge);

    M_tar = M_anat_amy;
    M_sm_amy= [];
    for ii = 1:length(histparam)-1
        for jj = 1:length(histparam)-1
            temp = M_tar(histparam(ii):histparam(ii+1),histparam(jj):histparam(jj+1));
            if ~isempty(temp)
                M_sm_amy(ii,jj) = mean(temp(:));
            end
        end
    end

    % plot
    figure()
    subplot(1,2,1);
    imagesc(M_sm_amy);
    caxis([0.97 1.03])
    colorbar
    title('Amy')
    subplot(1,2,2);
    imagesc(M_sm_ofc);
    caxis([0.97 1.03])
    colorbar
    title('ofc')

    % Valence RDM
    idx2 = std(M_val_img,1)==0;
    M_red = M_val_img(~idx2,~idx2);
    [idx] = kmeans(M_red,5,'Distance','correlation');
    [kidx] = kmeans(behav_ratings_,4);
    [~,argsort] =sort(kidx);
    imagesc( M_val(argsort,argsort))
    colormap(swampsunset)
    idx_exp = 1:1:length(M_val_img);
    idx_exp_km = idx_exp(idx);
    idx_exp_nkm = idx_exp(idx2);
    idx_idx = [idx_exp_km, idx_exp_nkm];
    M_val_img = M_val_img(idx_idx, idx_idx);
    imagesc(M_val_img)
end

%% RDMs
figure()
hold on
subplot(2,4,1)
plot(val_sc,val_sc)
ylabel('Valence')

subplot(2,4,2)
plot(val_sc,sal_sc)
ylabel('Salience')

subplot(2,4,3)
plot(val_sc,valp)
ylabel('Val+')

subplot(2,4,4)
plot(val_sc,valn)
ylabel('Val-')

subplot(2,4,5)
val_mat = 1-abs(val_sc-val_sc');
valp_mat = 1-abs(valp-valp');
valn_mat = 1-abs(valn-valn');
sal_mat = 1-abs(sal_sc-sal_sc');
imagesc(val_sc,val_sc,flipud(val_mat))
colormap(ARC_customcmap(64))
yticks([-1 -0.5 0 0.5 1])
yticklabels(-[-1 -0.5 0 0.5 1])
xlabel('Pleasantness')
ylabel('Pleasantness')

subplot(2,4,6)
imagesc(val_sc,val_sc,flipud(sal_mat))
colormap(ARC_customcmap(64))
yticks([-1 -0.5 0 0.5 1])
yticklabels(-[-1 -0.5 0 0.5 1])
xlabel('Pleasantness')

subplot(2,4,7)
imagesc(val_sc,val_sc,flipud(valp_mat))
colormap(ARC_customcmap(64))
yticks([-1 -0.5 0 0.5 1])
yticklabels(-[-1 -0.5 0 0.5 1])
xlabel('Pleasantness')

subplot(2,4,8)
imagesc(val_sc,val_sc,flipud(valn_mat))
colormap(ARC_customcmap(64))
yticks([-1 -0.5 0 0.5 1])
yticklabels(-[-1 -0.5 0 0.5 1])
xlabel('Pleasantness')
colorbar

%% P values of regression:

