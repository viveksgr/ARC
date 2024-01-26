%% General Settings
root = 'C:\Work\ARC\ARC';
maskfile =  'ARC3_anatgw.nii';
fmaskfile = 'ARC3_fanatgw3_pos.nii';
fmasker = false;
binz = 10;

anat_names = {'PC','AMY','OFC','VMPFC'};
anat_masks = {'rwPC.nii','rwAmygdala.nii','rwofc.nii','rwvmpfc.nii'};
% anat_names = {'PC','AMY','Thal','frontal_sal','frontal_val'}; 
% anat_masks = {'rwPC.nii','rwAmygdala.nii','rwThal.nii','frontal_sal.nii','frontal_val.nii'};
nanat = length(anat_names);
medianize_behav = true;
salier = true;
zscorer = true;
rangenormer = false;
sz_cntrl = false;

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
raw_RSA = false;     

nodor = 160;
intens_reg = false;
val_sep = true;
spearmanner = true;
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

if raw_RSA ; rsa_P1 = zeros(3,nanat,2); else; rsa_P1 = zeros(3,nanat,2); end % Complete

rs = zeros(3,1);
hold on
% Subject - index
kk = 1;
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
        
        % if zscorer
        %     modelmd2 = zscore(modelmd2,[],2);
        % end
        % [M1_new, M2_new] = ARC_transformMatrix(behav_ratings_,  modelmd2, binz,[ms1 ms2]);


        if ~rangenormer
            if ~sz_cntrl
                modelmd_binned = ARC_binAndTransform(modelmd2, behav_ratings_, 2*binz, [ms1 ms2]);
            else
                modelmd_binned = ARC_binAndTransform_sz(modelmd2, behav_ratings_, 2*binz, [ms1 ms2],5);
            end
        else
            modelmd_binned = ARC_binAndTransformQuantiles(modelmd2, behav_ratings_, 2*binz);
            edges_ = quantile(behav_ratings_, linspace(0, 1, 2*binz+1));
            centers_ = edges_(1:end-1)+(edges_(2)-edges_(1))/2;
            centers_rn = 2*((centers_-ms1)./(ms2-ms1))-1;
            centers_rd = round(centers_rn,2);
        end
        modelmd_corrcoef = corrcoef(modelmd_binned);
        
          
        M1_new = modelmd_binned(:,1:binz);
        M2_new = modelmd_binned(:,binz+1:end);         
        M1_corr = corrcoef(M1_new);
        M2_corr = corrcoef(M2_new);

        M_mat = corrcoef_2(M1_new',M2_new');


        % Make figure
        subplot(3,4,kk)
        hold on

        x = tsne(modelmd_binned');
        colormap = parula(20);
        % Create the scatter plot
        hold on;
        for i = 1:length(x)
            scatter(x(i,1), x(i,2), 36, colormap(i, :), 'filled');
        end

        % imagesc(M_mat)


        % if ~rangenormer
        %     xticks(0:5:20)
        %     yticks(0:5:20)
        %     yticklabels(-1:0.5:0)
        %     xticklabels(0:0.5:1)
        % else
        %     xticks(1:binz)
        %     yticks(1:binz)
        %     yticklabels(centers_rd(1:binz))
        %     xticklabels(centers_rd(binz+1:end))
        % end

        ylabel('Val-')
        xlabel('Val+')
        kk = kk+1;
        axis tight
        % colorbar
        % caxis([0.7 1 ])
        title(anat_names{ii})

        if raw_RSA
            utl_mask = logical(triu(ones(length( M2_corr)),1)); % All possible odors
            temp = corrcoef(M1_corr(utl_mask),M2_corr(utl_mask),'Rows','complete');
            % rsa_P1(s,ii,1) = temp(2);

            % antidiag = M_mat(sqrt(end):sqrt(end)-1:end-1);
            diager = diag(M_mat);

            % if salier
            %    M2_new = fliplr(M2_new);
            % end

            rsa_P1(s,ii,1) = nanmean(diager);
            M_mat = fliplr(M_mat);
            diager = diag(M_mat);
            rsa_P1(s,ii,2) = nanmean(diager);

            nel = M1_corr(utl_mask)+M2_corr(utl_mask);
            nel = sum(~isnan(nel));
            rs(s) = r2t(0.05,nel);

        else
            % if and(s==2,ii==2)
            %     'dil ki surkh deewaron pe'
            % end
            val_sc = linspace(-1,1,2*binz);
            sal_sc = abs(val_sc);
            valp = val_sc;
            valp(1:binz)=0;
            valn = val_sc;
            valn(binz+1:end)=0;
            

            val_mat = 1-abs(val_sc-val_sc');
            valp_mat = 1-abs(valp-valp');
            valn_mat = 1-abs(valn-valn');
            sal_mat = 1-abs(sal_sc-sal_sc');

            % [~, t_scores] = SFP_computeWeights_tsc_4reg( val_mat(:), sal_mat(:),  valp_mat(:),valn_mat(:),modelmd_corrcoef(:));
            % [~, t_scores] = SFP_computeWeights_tsc( valp_mat(:),valn_mat(:),modelmd_corrcoef(:));
            % [~, t_scores] = SFP_computeWeights_tsc( val_mat,sal_mat,modelmd_corrcoef,true);
            [~, t_scores] = SFP_computeWeights_tsc( val_mat,valn_mat,modelmd_corrcoef,true);
            % [~, t_scores] = SFP_computeWeights_tsc_3reg( val_mat,valp_mat,valn_mat,modelmd_corrcoef,true);
            % [~, t_scores] = SFP_computeWeights_tsc_4reg( val_mat,sal_mat,valp_mat,valn_mat,modelmd_corrcoef);
            rsa_P1(s,ii,:) = t_scores(2:end);

   
        end
    end
end
mkdir(savepath)

S_mat = squeeze(mean(rsa_P1));
S_err = squeeze(std(rsa_P1))./sqrt(3);

figure('Position',[0.5 0.5 640 480])
hold on
ngroups = size(S_mat, 1);
nbars = size(S_mat, 2);
b = bar(S_mat);
% b(1).FaceColor = [0 0.2470 0.9410];
% b(2).FaceColor = [0.3010 0.7450 0.9330];

% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
x_m = [];
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1)*groupwidth/(2*nbars);
    errorbar(x, S_mat(:,i), S_err(:,i), 'k.');
    x_m = [x_m; x];
end
xticks(1:nanat)
xticklabels(anat_names);
c_s = {'r','g','b'}; % Data dots for subjects
for ii = 1:nanat % For bars for perceptual, chemical and combinations
    for jj = 1:3
        plot(x_m(:,ii),squeeze(rsa_P1(jj,ii,:)),c_s{jj})
    end
end

clear modelmd_ modelmd modelmd2 S1_omat_vals S2_omat_vals unity M_anat M_sal M_val 
save(fullfile(savepath,'ARC_RSA'))

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