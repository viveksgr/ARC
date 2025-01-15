%% General Settings
root = 'C:\Work\ARC\ARC';
maskfile =  'ARC3_anatgw.nii';
fmaskfile = 'ARC3_fanatgw3_pos.nii';
fmasker = true;
binz =2 ;
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

dirs2 = {fullfile(root,'\ARC01\single');
    fullfile(root,'\ARC02\single');
    fullfile(root,'\ARC03\single')};

behav = load(fullfile(root,'ARC','NEMO_perceptual2.mat'));
% modelname = fullfile('mainval',sprintf('bin%01d',binz));
modelname = 'Identity_coding';
savepath = fullfile(root,'RSA',modelname);
v_id = 2; % Index of vector for median splitting odors

% load(fullfile(statpath,'fir_cv.mat'))
fprintf('\n')

rsa_P1t = zeros(3,nanat,2*binz);
rsa_P1t2 = zeros(3,nanat,2*binz);
rsa_P1wt = zeros(3,nanat,2*binz);
rsa_P1wt2 = zeros(3,nanat,2*binz);
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
load('C:\Work\ARC\ARC\RSA\main_perm\ARC_RSA.mat','t_score_mat')

for s = [1 2 3] % Subject
    fprintf('Subject: %02d\n',s)
    anatdir = fullfile(root,sprintf('ARC%02d',s),'single');
    % Construction of median split
    behav_ratings = behav.behav(s).ratings(:,v_id);
    behavmat =  behav.behav(s).ratings(:,setdiff(1:18,v_id));
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

    path2 = dirs2{s};
    load(fullfile(path2 ,'task_struct_trialwise.mat'))
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
    behav_matrix = corrcoef(behavmat(group_vec,:)');
    behav_ratings_disc = discretize(behav_ratings_ ,binz);

    for ii = 1:length(anat_names)
        fprintf('area:%02d\n',ii)
        modelmd_ = load(fullfile(anatdir,anat_names{ii},'TYPED_FITHRF_GLMDENOISE_RR.mat'),'modelmd','noisepool');
        modelmd = squeeze(modelmd_.modelmd);
        noisepool = modelmd_.noisepool;
        if single_c
            modelmd = modelmd(masks_set_cell{ii},:);
        end
        [r1,~] = find(isnan(modelmd));
        modelmd(r1,:) = [];
        modelmd = zscore(modelmd,[],2);

        % Binwise RSA for whole population
        % anat_corr = corrcoef(modelmd);
        % for pp = 1:binz
        %     trial_disc = behav_ratings_disc==pp;
        %     trial_mat = and(trial_disc,trial_disc');
        %     trial_mask = and(trial_mat,utl_mask);
        % 
        %     [temp,tsc] =  ARC_multicomputeWeights_tsc( [behav_matrix(trial_mask) unity(trial_mask) task_run(trial_mask)...
        %                   sess_run(trial_mask) set_run(trial_mask)], anat_corr(trial_mask)) ;
        %     rsa_P1wt(s,ii,pp) = temp(2);
        %     rsa_P1wt2(s,ii,pp) = temp(3);
        %     rsa_P1t(s,ii,pp) = tsc(2);
        %     rsa_P1t2(s,ii,pp) = tsc(3);
        %     [rsa_Pcorr(s,ii),t_scores, rsa_prop(s,ii,:),w_scores] = ARC_multicomputeWeights_tsc_voxwise(valp_mat, valn_mat, modelmd_binned);
        % end

        assert(size(modelmd,1)==size(t_score_mat{s,ii},1))
        thr = tinv(0.975,size(t_score_mat{s,ii},1));
        pospop = and(t_score_mat{s,ii}(:,1)>thr,t_score_mat{s,ii}(:,2)<thr);
        negpop = and(t_score_mat{s,ii}(:,2)>thr,t_score_mat{s,ii}(:,1)<thr);
        % pospop = t_score_mat{s,ii}(:,1)>thr;
        % negpop = t_score_mat{s,ii}(:,2)>thr;
        

        vid_idx = behav_ratings_>median(behav_ratings_);
        vid_mat_pos = and(vid_idx,vid_idx');
        vid_mat_neg = and(~vid_idx,~vid_idx');


        trial_mask = and(utl_mask,vid_mat_pos);
        anat_corr = corrcoef(modelmd(pospop,:));
        [temp,tsc] =  ARC_multicomputeWeights_tsc( [behav_matrix(trial_mask) unity(trial_mask) task_run(trial_mask)...
                          sess_run(trial_mask) set_run(trial_mask)], anat_corr(trial_mask)) ;
        rsa_P1wt(s,ii,1) = temp(2);
        rsa_P1wt2(s,ii,1) = temp(3);
        rsa_P1t(s,ii,1) = tsc(2);
        rsa_P1t2(s,ii,1) = tsc(3);
        
        anat_corr = corrcoef(modelmd(negpop,:));
        [temp,tsc] =  ARC_multicomputeWeights_tsc( [behav_matrix(trial_mask) unity(trial_mask) task_run(trial_mask)...
                          sess_run(trial_mask) set_run(trial_mask)], anat_corr(trial_mask)) ;
        rsa_P1wt(s,ii,2) = temp(2);
        rsa_P1wt2(s,ii,2) = temp(3);
        rsa_P1t(s,ii,2) = tsc(2);
        rsa_P1t2(s,ii,2) = tsc(3);



        trial_mask = and(utl_mask,vid_mat_neg);
        anat_corr = corrcoef(modelmd(pospop,:));
        [temp,tsc] =  ARC_multicomputeWeights_tsc( [behav_matrix(trial_mask) unity(trial_mask) task_run(trial_mask)...
                          sess_run(trial_mask) set_run(trial_mask)], anat_corr(trial_mask)) ;
        rsa_P1wt(s,ii,3) = temp(2);
        rsa_P1wt2(s,ii,3) = temp(3);
        rsa_P1t(s,ii,3) = tsc(2);
        rsa_P1t2(s,ii,3) = tsc(3);
        
        anat_corr = corrcoef(modelmd(negpop,:));
        [temp,tsc] =  ARC_multicomputeWeights_tsc( [behav_matrix(trial_mask) unity(trial_mask) task_run(trial_mask)...
                          sess_run(trial_mask) set_run(trial_mask)], anat_corr(trial_mask)) ;
        rsa_P1wt(s,ii,4) = temp(2);
        rsa_P1wt2(s,ii,4) = temp(3);
        rsa_P1t(s,ii,4) = tsc(2);
        rsa_P1t2(s,ii,4) = tsc(3);







    end
end
p_values_3d = ARC_RSA_pvals(rsa_P1t2, rsa_P1wt2, sum(trial_mask(:)))
p_values_matrix = ARC_RSA_pvals_diff(cat(3,rsa_P1t,rsa_P1t2), cat(3,rsa_P1wt,rsa_P1wt2), sum(trial_mask(:)))

ARC_barplot(rsa_P1wt2)
ylabel('RSA wt')
xticks(1:4)
xticklabels(anat_names)
mkdir(savepath)
legend({'Val+ in +ve vox','Val+ in -ve vox','Val- in +ve vox','Val- in -ve vox'})
savefig(fullfile(savepath,'ARC_RSAwt_ident'))
print(fullfile(savepath,'ARC_RSAwt_ident'),'-dpng')

clear modelmd_ modelmd modelmd2 S1_omat_vals S2_omat_vals unity M_anat M_sal M_val
save(fullfile(savepath,'ARC_RSA'))


