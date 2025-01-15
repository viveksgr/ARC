%% General Settings
root = 'C:\Work\ARC\ARC';
maskfile =  'ARC3_anatgw.nii';
fmaskfile = 'ARC3_fanatgw3_pos.nii';
fmasker = false; % Turn off the functional mask
anat_names = {'PC','AMY','OFC','vmpfc'};
anat_masks = {'rwPC.nii','rwAmygdala.nii','rwofc.nii','rwvmpfc.nii'};
nanat = length(anat_names);
intens_reg = false;
dirs = {fullfile(root,'\ARC01\mediation');
    fullfile(root,'\ARC02\mediation');
    fullfile(root,'\ARC03\mediation')};

behavP = load(fullfile(root,'ARC','NEMO_perceptual2.mat'));
behavC = load(fullfile(root,'ARC','NEMO_chemical2.mat'));

modelname = 'searchl_standard';
loadpath = fullfile('C:\Work','\ARC\ARC\RSA',modelname); % Where anat files are stored
modelname2 = 'searchl_standard\sal';
savepath = fullfile('C:\Work','\ARC\ARC\RSA',modelname2);
mkdir(savepath)
v_id = 2; % Index of vector for median splitting odors

valier = false;
main_masker = true;
featurer = false;
singletrialer = true;
percepter = true; % Use perceptual rather than chemical ratings

nodor = 160;
single_n = false; % Noisepool
single_c = true; % Cutoff from sign voxels

% load(fullfile(statpath,'fir_cv.mat'))
fprintf('\n')
rsa_P1_ = zeros(3,nchoosek(nanat,2)); % Complete
hold on
anat_cell = {};
% Subject - index
figure('Position',[0.5 0.5 1280 320])
for s = [1 2 3] % Subject
    fprintf('Subject: %02d\n',s)
    anatdir = fullfile(root,sprintf('ARC%02d',s),'single');
    % Construction of median split
    behav_ratingsC = behavC.behav(s).ratings(:,1:20); % take only first 20 components of chemical space
    behav_ratingsP = behavP.behav(s).ratings(:,2); % 
    if percepter; X_mat = behavP.behav(s).ratings(:,[1 3:end]); else; X_mat =  behav_ratingsC; end

    [~, ~, p_values_val] = bootstrapRidge(X_mat, behav_ratingsP, 1000, 0.1);
    [~, ~, p_values_sal] = bootstrapRidge(X_mat, abs(behav_ratingsP), 1000, 0.1);

    %% Gray Matter, Functional and Anatomical Masks
    statpath = dirs{s};
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
    anatpath = anatdir;
        
    % Model names
    masks_set = [];
    masks_set_cell = {};
    
    anatmasks = [];
    for ii = 1:length(anat_masks)
        m1 = spm_read_vols(spm_vol(fullfile(anatpath,anat_masks{ii})));
        m1(isnan(m1))=0;
        m1(m1<=0.01)=0;
        m1(m1>0) = 1;
        m1 = m1(mask);
        masks_set(:,ii)=m1(fmask_1d);
        
        masks_set_cell{ii} = fmask_1d(logical(m1));
        anatmask = (spm_read_vols(spm_vol(fullfile(anatpath, anat_masks{ii}))));
        anatmask(isnan(anatmask)) = 0;
        anatmask = logical(anatmask);
        anatmask = and(anatmask,marea);
        anatmasks(:,:,:,ii) = anatmask;
    end
    anat_cell{s}= anatmasks;
    masks_set(isnan(masks_set))=0;
    %%{
    %% Defining RDMs
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
    n_reg = unity(logical(utl_mask));
    if singletrialer; Pmat_group = X_mat(group_vec,:) ;  else Pmat_group = X_mat; end

    if main_masker; thr = r2t(0.05,nchoosek(length(group_vec),2)); else; thr = -99; end

    if valier
        Pmat_val = Pmat_group(:,p_values_val<0.05);
        behavioral_corr = corrcoef(Pmat_val');
        behavioral_corr_lin = behavioral_corr(utl_mask);
        behavioral_corr_lin(isnan(behavioral_corr_lin))=1;
        valcorr = -abs(behav_ratingsP(group_vec)-behav_ratingsP(group_vec)');
        valcorr_lin = valcorr(utl_mask);
    else
        Pmat_sal = Pmat_group(:,p_values_sal<0.05);
        behavioral_corr = corrcoef(Pmat_sal');
        behavioral_corr_lin = behavioral_corr(utl_mask);
        behavioral_corr_lin(isnan(behavioral_corr_lin))=1;
        valcorr = -abs(abs(behav_ratingsP(group_vec))-abs(behav_ratingsP(group_vec))');
        valcorr_lin = valcorr(utl_mask);
    end

    if featurer; regvec = behavioral_corr_lin; fpre = 'P_ARC'; else; regvec = valcorr_lin; fpre = 'ARC';end
    %% Representational connectivity
    anat_corr_mat = {};
    for ii = 1:length(anat_names)
        fprintf('area:%02d\n',ii)
        modelmd_ = load(fullfile(anatdir,anat_names{ii},'TYPED_FITHRF_GLMDENOISE_RR.mat'),'modelmd','noisepool');
        modelmd = squeeze(modelmd_.modelmd);
        noisepool = modelmd_.noisepool;
        if single_c
            modelmd = modelmd(masks_set_cell{ii},:);
            noisepool = noisepool(masks_set_cell{ii});
        end
        S_omat_vals = zscore(modelmd,[],2);
        [r1,~] = find(isnan(S_omat_vals));
        S_omat_vals(r1,:) = [];
    
        % fileload for finding voxels with significant performance
        if valier
            fmat = spm_read_vols(spm_vol(sprintf('%s%02d_val_RSAx1%s.nii',fpre,s,anat_names{ii})));
        else
            fmat = spm_read_vols(spm_vol(sprintf('%s%02d_sal_RSAx1%s.nii',fpre,s,anat_names{ii})));
        end        
        fmat_mask = fmat(logical(anatmasks(:,:,:,ii)));
        assert(size(fmat_mask,1)==size(S_omat_vals,1),'Mask and Anat areas differ in size')

        fmat_mask_thr = fmat_mask>thr;
        anat_corr_mat{ii} = corrcoef(S_omat_vals(fmat_mask_thr,:));
    end

    %% RCA
    kk = 0;
    for ii = 1:nanat
        for jj = ii+1:nanat
            kk = kk+1;
            a1 = anat_corr_mat{ii};
            a2 = anat_corr_mat{jj};
            a1 = a1(utl_mask);
            a2 = a2(utl_mask);
            a1 = regressmeout(a1',n_reg')';
            a2 = regressmeout(a2',n_reg')';
            
            a1 = a1.* regvec ;
            a2 = a2.* regvec; 
            tmp = corrcoef(a1,a2);
            % tmp2 = partialcorr(a1,a2, valcorr_lin );

             rsa_P1_(s,kk) = tmp(2);
        end
    end

    % subplot(1,3,s)
    % imagesc(rsa_mat_subject)%,'alphadata',0.3+0.7*double(abs(rsa_mat_subject)>thr))
    % xticks(1:length(anat_names))
    % xticklabels(anat_names)
    % xtickangle(90)
    % yticks(1:length(anat_names))
    % yticklabels(anat_names)
    % colorbar
%     caxis([-0.2 0.2])
    if valier; title(sprintf('Subject: %02d Valence',s)); else title(sprintf('Subject: %02d Salience',s)); end  

end
save(fullfile(savepath,"rsa_P1"),'rsa_P1_')
savefig(fullfile(savepath,'ARC_RSA_conn_sub.fig'))
print(fullfile(savepath,'ARC_RSA_conn_sub'),'-dpng')

% figure()
% mmat = squeeze(mean(rsa_P1_,1));
% imagesc(mmat)%,'alphadata',0.3+0.7*double(abs(rsa_mat_subject)>thr))
% xticks(1:length(anat_names))
% xticklabels(anat_names)
% xtickangle(90)
% yticks(1:length(anat_names))
% yticklabels(anat_names)
% colorbar
% %     caxis([-0.2 0.2])
% if valier; title('Val'); else title('Sal'); end
% savefig(fullfile(savepath,'ARC_RSA_conn_grp.fig'))
% print(fullfile(savepath,'ARC_RSA_conn_grp'),'-dpng')

%% Bar plots
anat_names2 = {'PC-AMY','PC-OFC','PC-VMPFC','AMY-OFC','AMY-VMPFC','OFC-VMPFC'};
dirs = {'C:\Work\ARC\ARC\RSA\searchl_standard\val';'C:\Work\ARC\ARC\RSA\searchl_standard\sal'};
rsa_P1 = variable_extract(dirs,'rsa_P1.mat','rsa_P1_',false);
rsa_P1 = cat(3,rsa_P1{:});
multisubbarplot(rsa_P1,anat_names2,'RCA')
legend({'Salience','Salience features'},'Location','northwest')
ylim([0 0.35])
savefig("RCA_sal")
print("RCA_sal","-dpng")


