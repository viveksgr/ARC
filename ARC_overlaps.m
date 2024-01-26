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
loadpath1 = fullfile('C:\Work','\ARC\ARC\RSA',modelname); % Where anat files are stored

modelname = 'searchl_perc';
loadpath2 = fullfile('C:\Work','\ARC\ARC\RSA',modelname); % Where anat files are stored

modelname2 = 'RCA';
savepath = fullfile('C:\Work','\ARC\ARC\RSA',modelname2);
mkdir(savepath)
v_id = 2; % Index of vector for median splitting odors

singletrialer = true;

nodor = 160;
single_n = false; % Noisepool
single_c = true; % Cutoff from sign voxels

% load(fullfile(statpath,'fir_cv.mat'))
fprintf('\n')
rsa_P1_ = zeros(3,nanat,2); % Complete
hold on
anat_cell = {};
% Subject - index
thr = 5.3853e-04;
for s = [1 2 3] % Subject
    fprintf('Subject: %02d\n',s)
    anatdir = fullfile(root,sprintf('ARC%02d',s),'single');
    % Construction of median split
 
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

        fmat_val = spm_read_vols(spm_vol(fullfile(loadpath1,sprintf('ARC%02d_val_RSAx1%s.nii',s,anat_names{ii}))));
        fmat_val_feat = spm_read_vols(spm_vol(fullfile(loadpath2,sprintf('P_ARC%02d_val_RSAx1%s.nii',s,anat_names{ii}))));

        fmat_sal = spm_read_vols(spm_vol(fullfile(loadpath1,sprintf('ARC%02d_sal_RSAx1%s.nii',s,anat_names{ii}))));
        fmat_sal_feat = spm_read_vols(spm_vol(fullfile(loadpath2,sprintf('P_ARC%02d_sal_RSAx1%s.nii',s,anat_names{ii}))));

        fmat_val = fmat_val(logical(anatmasks(:,:,:,ii)))>thr;
        fmat_val_feat = fmat_val_feat(logical(anatmasks(:,:,:,ii)))>thr;
        fmat_sal = fmat_sal(logical(anatmasks(:,:,:,ii)))>thr;
        fmat_sal_feat = fmat_sal_feat(logical(anatmasks(:,:,:,ii)))>thr;

        rsa_P1_(s,ii,1) = sum(and(fmat_val,fmat_val_feat))/sum(fmat_val);
         rsa_P1_(s,ii,2) = sum(and(fmat_sal,fmat_sal_feat))/sum(fmat_sal);

    end
  

end
save(fullfile(savepath,"rsa_P1"),'rsa_P1_')


multisubbarplot(  rsa_P1_,anat_names,'Overlap')
legend({'Valence','Salience'},'Location','northwest')
% ylim([0 0.35])
savefig("RCA_overlap_stand")
print("RCA_overlap_stand","-dpng")


