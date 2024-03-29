%% General Settings
root = 'C:\Work\ARC\ARC';
maskfile =  'ARC3_anatgw.nii';
fmaskfile = 'ARC3_fanatgw3_pos.nii';
fmasker = false;
binz = 3;
if mod(binz,2)==0; binzpart1 = binz/2; binzpart2 = binzpart1+1; else; binzpart1 = (binz+1)/2 ; binzpart2 = binzpart1; end

anat_names = {'PC','AMY','OFC','VMPFC'};
anat_masks = {'rwPC.nii','rwAmygdala.nii','rwofc.nii','rwvmpfc.nii'};
% anat_names = {'PC','AMY','Thal','frontal_sal','frontal_val'}; 
% anat_masks = {'rwPC.nii','rwAmygdala.nii','rwThal.nii','frontal_sal.nii','frontal_val.nii'};
nanat = length(anat_names);
medianize_behav = true;
salier = true;
zscorer = true;
rangenormer = true;
sz_cntrl = false;
raw_RSA = false;     

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

spearmanner = true;
sess_l = repmat([0],1,2,3);
dirs = {fullfile(root,'\ARC01\mediation');
    fullfile(root,'\ARC01\mediation');
    fullfile(root,'\ARC01\mediation')};
behav = load(fullfile(root,'ARC','NEMO_perceptual2.mat'));
modelname = 'Decoding';
savepath = fullfile(root,'RSA',modelname);
v_id = 2; % Index of vector for median splitting odors

% load(fullfile(statpath,'fir_cv.mat'))
fprintf('\n')

rsa_P1 = zeros(3,nanat,2); 
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
        
        if rangenormer
            edges_val = quantile(behav_ratings_,linspace(0,1,binz+1));
            edges_sal = quantile(abs(behav_ratings_),linspace(0,1,binz+1));
            labels_val = discretize(behav_ratings_,edges_val);
            labels_sal = discretize(abs(behav_ratings_),edges_sal);
        else
            labels_val = discretize(behav_ratings_,binz);
            labels_sal = discretize(abs(behav_ratings_),binz);
        end

        rsa_P1(s,ii,1) = Classify_Permute_VS2( double(modelmd2'), labels_val, 10);
        rsa_P1(s,ii,2) = Classify_Permute_VS2( double(modelmd2'), labels_sal, 10);

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
ylabel('Model accuracy')
legend({'Val','Sal'})
yline(1/binz)
% legend({'Valence','Salience'})
clear modelmd_ modelmd modelmd2 S1_omat_vals S2_omat_vals unity M_anat M_sal M_val 

p_values = arrayfun(@(r_eff) ARC_computePValueOneTailed(r_eff, binz, length(behav_ratings_)), rsa_P1);

savefig(fullfile(savepath,'ARC_decoding'))
save(fullfile(savepath,'ARC_decoding'))