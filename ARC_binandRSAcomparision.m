%% General Settings
root = 'C:\Work\ARC\ARC';
maskfile =  'ARC3_anatgw.nii';
fmaskfile = 'ARC3_fanatgw3_pos.nii';
fmasker = true;

anat_names = {'PC','AMY','OFC','VMPFC'};
anat_masks = {'rwPC.nii','rwAmygdala.nii','rwofc.nii','rwvmpfc.nii'};
% anat_names = {'PC','AMY','Thal','frontal_sal','frontal_val'};
% anat_masks = {'rwPC.nii','rwAmygdala.nii','rwThal.nii','frontal_sal.nii','frontal_val.nii'};
nanat = length(anat_names);
medianize_behav = false;
% sess_l = cat(3,nchoosek([1 2 3],2),nchoosek([2 3 4],2),nchoosek([2 3
% 4],2),nchoosek([2 3 4],2)); % For sesswise
% load('C:\Data\NEMO\swampsunset.mat');

single_trial = true;
single_z = true; % Zscore
coarser_ = false; % Only binary estimation of valence RSM
single_n = false; % Noisepool
single_c = true; % Cutoff from sign voxels
increm = false; % Compare rating to previous
partial_corr = false;

nodor = 160;
intens_reg = false;
val_sep = false;
spearmanner = false;
sess_l = repmat([0],1,2,3);
dirs = {fullfile(root,'\ARC01\mediation');
    fullfile(root,'\ARC01\mediation');
    fullfile(root,'\ARC01\mediation')};
behav = load(fullfile(root,'ARC','NEMO_perceptual2.mat'));
modelname = 'temp';
savepath = fullfile(root,'RSA',modelname);
v_id = 2; % Index of vector for median splitting odors

% load(fullfile(statpath,'fir_cv.mat'))
fprintf('\n')
rsa_P1 = zeros(3,nanat,2); % Complete
rsa_P2 = zeros(3,nanat,2);
hold on
% Subject - index
for s = [1 2 3] % Subject
    fprintf('Subject: %02d\n',s)
    anatdir = fullfile(root,sprintf('ARC%02d',s),'single');
    % Construction of median split
    behav_ratings = behav.behav(s).ratings(:,v_id);
    behav_int = behav.behav(s).ratings(:,1);
    % behav_ratings = vs_normalizer(behav_ratings);
    if medianize_behav
        behav_ratings = normalize(behav_ratings,'medianiqr');
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
    if s <3
        anatpath = fullfile(root,sprintf('NEMO_%02d',s),'\imaging\nii\masks');
    else
        anatpath = fullfile(root,sprintf('NEMO_%02d',s+1),'\imaging\nii\masks');
    end

    % Model names
    masks_set = [];
    masks_set_cell = {};
    for ii = 1:length(anat_masks)
        m1 = spm_read_vols(spm_vol(fullfile(anatdir,anat_masks{ii})));
        m1(isnan(m1))=0;
        m1(m1<=0.01)=0;
        m1(m1>0) = 1;
        m1 = m1(mask);
        masks_set(:,ii)=m1(fmask_1d);
        fprintf('area count: %04d\n',sum(m1(fmask_1d)))
        if single_trial
            masks_set_cell{ii} = fmask_1d(logical(m1));
        end
    end
    masks_set(isnan(masks_set))=0;
    linux_config = false;
    warning('off','all')

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
    bv = abs(behav_ratings_);
    if coarser_
        behav_ratings_(behav_ratings_<=md) = -1;
        behav_ratings_(behav_ratings_>md) = 1;
        md2 = median(bv);
        bv(bv<=md2) = -1;
        bv(bv>md2) = 1;
    end

    M_val = abs(behav_ratings_-behav_ratings_');
    %     M_val = -(behav_ratings_.*behav_ratings_');
    [~,argsort] = sort(behav_ratings_);

    %         figure('Position',[0.5 0.5 300 240])
    %         M_val_img = M_val(argsort,argsort);
    %         imagesc(M_val_img)
    %         colorbar
    %         caxis([0 3])
    %         title('Valence')
    %         colormap(swampsunset)
    %         print(fullfile(savepath,sprintf('Val%02d',s)),'-dsvg')


    M_val = M_val(utl_mask);
    M_sal = abs(bv-bv');
    %         figure('Position',[0.5 0.5 300 240])
    %         M_val_img = M_sal(argsort2,argsort2);
    %         imagesc(M_val_img)
    %         colorbar
    %         caxis([0 1.5])
    %         title('Salience')
    %         colormap(swampsunset)
    %         print(fullfile(savepath,sprintf('Sal%02d',s)),'-dsvg')


    M_sal = M_sal(utl_mask);


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


        if ii==3
            'loopdeloop'
        end

        if and(single_z,~val_sep)
            S1_omat_vals = zscore(modelmd2,[],2);
        else
            S1_omat_vals = modelmd2;
        end


        [r1,~] = find(isnan(S1_omat_vals));
        S1_omat_vals(r1,:) = [];
        M_anat = 1-corrcoef(S1_omat_vals);
        %         M_anat = pdist(S1_omat_vals','Mahalanobis');
        %         M_anat = squareform(M_anat);



        if intens_reg
            behav_int =  behav_int(group_vec);
            M_int = -( behav_int.* behav_int');
            M_int = M_int(utl_mask);
            M_val = regressmeout(M_val', M_int')';
            M_sal = regressmeout(M_sal', M_int')';
        end

        M_anat = M_anat(utl_mask);
        if partial_corr
            M_anat = regressmeout(M_anat', unity(utl_mask)')';
            %                 temp = partialcorr(M_anat(utl_mask),M_val,unity(utl_mask));
            %                 temp2 = partialcorr([M_anat(utl_mask),M_sal],unity(utl_mask));
            %             else
        end

        temp = corrcoef(M_anat,M_val);
        temp2 = corrcoef(M_anat,M_sal);
        if spearmanner
            temp = corr(M_anat,M_val,'type','spearman');
            temp(2) = temp;
            temp2 = corr(M_anat,M_sal,'type','spearman');
            temp2(2) = temp2;
        end

        S_mat(ii,1) = temp(2);
        [S_mat2(ii,1)] = r2p(double(temp(2)),length(M_val));
        S_mat(ii,2) = temp2(2);
        [S_mat2(ii,2)] = r2p(double(temp2(2)),length(M_sal));
        utl_mask2 = utl_mask;


    end
    rsa_P1(s,:,:) = S_mat;
    rsa_P2(s,:,:) = S_mat2;
end

%% Testing and debugging binning vs unbinned RSA

% neuralact and behav_ratings_ are given

% --- Unbinned-case
neuralact_z = zscore(neuralact,[],2);
corr_neural = corrcoef(neuralact_z);

% Compute valence and salience
salience_ratings = abs(behav_ratings_);
corr_salience = -abs(salience_ratings-salience_ratings');
corr_valence = -abs(behav_ratings_-behav_ratings_');

% Extract upper triangle
utl_mask = logical(triu(ones(length(corr_salience)),1)); % All possible odors

% raw_RSA = corrcoef(corr_neural(utl_mask),corr_valence(utl_mask));
% [~,t_st] = r2p(raw_RSA(2),sum(utl_mask(:)));

% Compute t-score in linear regression
[~, t_st1] = ARC_computeWeights_tsc(corr_valence(utl_mask), corr_salience(utl_mask),corr_neural(utl_mask));


% --- Binned-case
nbins = 60;

% Bin the neural activity
neuralact_binned = ARC_binAndTransform(neuralact, behav_ratings_, 2*nbins, [-1 1]);
corr_neural_binned = corrcoef(neuralact_binned);

% Create theoretical valence and salience
val_sc = linspace(-1,1,2*nbins);
sal_sc = abs(val_sc);
val_mat = 1-abs(val_sc-val_sc');
sal_mat = 1-abs(sal_sc-sal_sc');

% Extract upper triangle
utl_mask_binned = logical(triu(ones(length(sal_mat)),1)); % All possible odors

% Compute t-score in linear regression
[~, t_st2] = ARC_computeWeights_tsc(val_mat(utl_mask_binned), sal_mat(utl_mask_binned),corr_neural_binned(utl_mask_binned));


%% Testing and debugging binning vs unbinned RSA - 
% neuralact and behav_ratings_ are given

% --- Unbinned-case
neuralact_z = zscore(neuralact,[],2);
corr_neural = corrcoef(neuralact_z);

% Compute valence and salience
salience_ratings = abs(behav_ratings_);
corr_salience = -abs(salience_ratings-salience_ratings');
corr_valence = -abs(behav_ratings_-behav_ratings_');

% Extract upper triangle
utl_mask = logical(triu(ones(length(corr_salience)),1)); % All possible odors

corr_neural = corr_neural(utl_mask);
corr_neural = regressmeout(corr_neural', unity(utl_mask)')';

raw_RSA1 = corrcoef(corr_neural,corr_valence(utl_mask));
% [~,t_st] = r2p(raw_RSA(2),sum(utl_mask(:)));
raw_RSA2 = corrcoef(corr_neural,corr_salience(utl_mask));
% [~,t_st2] = r2p(raw_RSA(2),sum(utl_mask(:)));

% Non-linear
raw_RSA11 = corr(corr_neural,corr_valence(utl_mask),"type",'Spearman');
% [~,t_st_] = r2p(raw_RSA,sum(utl_mask(:)));
raw_RSA22 = corr(corr_neural,corr_salience(utl_mask),"type",'Spearman');
% [~,t_st2_] = r2p(raw_RSA,sum(utl_mask(:)));
bar([raw_RSA11-raw_RSA1(2) raw_RSA22-raw_RSA2(2)])

% Compute t-score in linear regression
[~, t_st3] = ARC_computeWeights_tsc(corr_valence(utl_mask), corr_salience(utl_mask),corr_neural);


% --- Binned-case
nbins = 10;

% Bin the neural activity
neuralact_binned = ARC_binAndTransform(neuralact, behav_ratings_, 2*nbins, [-1 1]);
corr_neural_binned = corrcoef(neuralact_binned);

% Create theoretical valence and salience
val_sc = linspace(-1,1,2*nbins);
sal_sc = abs(val_sc);
val_mat = 1-abs(val_sc-val_sc');
sal_mat = 1-abs(sal_sc-sal_sc');

% Extract upper triangle
utl_mask_binned = logical(triu(ones(length(sal_mat)),1)); % All possible odors

% Compute t-score in linear regression
[~, t_st4] = ARC_computeWeights_tsc(val_mat(utl_mask_binned), sal_mat(utl_mask_binned),corr_neural_binned(utl_mask_binned));

% Binning iteration
binz = 1:80;
t_st_mat = zeros(2,length(binz));
kk = 0;
for nbins = binz
    kk = kk+1;
    % Bin the neural activity
    neuralact_binned = ARC_binAndTransform(neuralact, behav_ratings_, 2*nbins, [-1 1]);
    corr_neural_binned = corrcoef(neuralact_binned);

    % Create theoretical valence and salience
    val_sc = linspace(-1,1,2*nbins);
    sal_sc = abs(val_sc);
    val_mat = 1-abs(val_sc-val_sc');
    sal_mat = 1-abs(sal_sc-sal_sc');

    % Extract upper triangle
    utl_mask_binned = logical(triu(ones(length(sal_mat)),1)); % All possible odors

    % Compute t-score in linear regression
    [~, t_st4] = ARC_computeWeights_tsc(val_mat(utl_mask_binned), sal_mat(utl_mask_binned),corr_neural_binned(utl_mask_binned));
    t_st_mat(:,kk) = t_st4(2:end);
end

plot(binz*2,t_st_mat')