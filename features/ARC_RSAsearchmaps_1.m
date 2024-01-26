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

modelname = 'standard_searchl';
savepath = fullfile(root,'\ARC\RSA',modelname);
mkdir(savepath)

v_id = 2; % Index of vector for median splitting odors
dister = false; % Use euclidean dist
singletrialer = false;
percepter = true; % Use perceptual rather than chemical ratings
intreg = true; % Regress out intensity
nodor = 160;
mapper  = true;

single_n = false; % Noisepool
single_c = true; % Cutoff from sign voxels

% load(fullfile(statpath,'fir_cv.mat'))
fprintf('\n')
rsa_P1_ = cell(3,nanat,2); % Complete
hold on
anat_cell = {};
% Subject - index
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
    if singletrialer; Pmat_group = X_mat(group_vec,:) ;  else Pmat_group = X_mat; end
    
    %% Representational connectivity
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
        %         M_anat = 1-corrcoef(S1_omat_vals);
        
        % Searchlight configuration
        rget = 3;
        sz = size(logical(anatmasks(:,:,:,ii)));
        ref_vox = round(sz/2);
        [MX, MY, MZ] = ndgrid(1:sz(1),1:sz(2),1:sz(3));
        radii = sqrt((MX-ref_vox(1)).^2 + (MY-ref_vox(2)).^2 + (MZ-ref_vox(3)).^2);
        % prototype sphere index for radii<rget that can be used everywhere in the brain
        radius_index = find(radii<rget) - sub2ind(sz,ref_vox(1),ref_vox(2),ref_vox(3)); % Neighbors index
        
        % Gray mask configuration
        lin_index = find(logical(anatmasks(:,:,:,ii)));
        nvox=length(lin_index);
        % resultsvolume index
        linindexconv = zeros(sz);
        linindexconv(lin_index) = 1:length(lin_index);
        
        rsa_vec_val = zeros(nvox,1);
        rsa_vec_sal = zeros(nvox,1);
        for cnt2 = 1:nvox
            indexindex2 = radius_index + lin_index(cnt2);
            indexindex2 = intersect(lin_index, indexindex2);
            indexindex2 = linindexconv(indexindex2);
            S_omat_vals_r = S_omat_vals(indexindex2,:);
            [r1,~] = find(isnan(S_omat_vals_r));
            S_omat_vals_r(r1,:) = [];
            if ~singletrialer;  S_omat_vals_r = splitapply( @mean,S_omat_vals_r',group_vec)'; end
            M_anat = 1-corrcoef(S_omat_vals_r);
        
            if (size(S_omat_vals_r,1)>1)
                Pmat_val = Pmat_group(:,p_values_val<0.05);
                Pmat_sal = Pmat_group(:,p_values_sal<0.05);

                if intreg
                    assert(percepter==true)
                end
                if sum(p_values_val<0.05)>=1
                    rsa_vec_val(cnt2,1) = correlateDistanceMatrix2(M_anat, Pmat_val,Pmat_group(:,1),intreg);
                end
                if sum(p_values_sal<0.05)>=1
                    rsa_vec_sal(cnt2,1) = correlateDistanceMatrix2(M_anat, Pmat_sal,Pmat_group(:,1),intreg);
                end
            end
        end
        rsa_P1_{s,ii,1} = rsa_vec_val;
        rsa_P1_{s,ii,2} = rsa_vec_sal;
        
        if mapper
            rsa_vec_val = unmasker(rsa_vec_val,logical(anatmasks(:,:,:,ii)));
            rsa_vec_sal = unmasker(rsa_vec_sal,logical(anatmasks(:,:,:,ii)));
            if percepter
                write_reshaped_nifty(rsa_vec_val, savepath, false, fullfile(statpath,'ARC3_anatgw.nii'), sprintf('P_ARC%02d_val_RSAx1%s',s,anat_names{ii}));
                write_reshaped_nifty(rsa_vec_sal, savepath, false, fullfile(statpath,'ARC3_anatgw.nii'), sprintf('P_ARC%02d_sal_RSAx1%s',s,anat_names{ii}));
            else
                write_reshaped_nifty(rsa_vec_val, savepath, false, fullfile(statpath,'ARC3_anatgw.nii'), sprintf('C_ARC%02d_val_RSAx1%s',s,anat_names{ii}));
                write_reshaped_nifty(rsa_vec_sal, savepath, false, fullfile(statpath,'ARC3_anatgw.nii'), sprintf('C_ARC%02d_sal_RSAx1%s',s,anat_names{ii}));
            end
        end
    end
    %}
end
clear M_val unity M_anat modelmd_ modelmd S_omat_vals utl_mask utl_mask2
save(fullfile(savepath,'ARC_RSA'))

ntrials = 160;
rsa_P1_pvals = cellfun(@(x) r2p(x,nchoosek(ntrials,2)),rsa_P1_,'UniformOutput',false);

rsa_P1_pvals2 = cat(1,rsa_P1_pvals{:});
rsa_P1_thresh = fdr_benjhoc(rsa_P1_pvals2);

% for pp = 1:numel(rsa_P1_thresh);  if isempty(rsa_P1_thresh{pp}); rsa_P1_thresh{pp} = 0; end; end

% rsa_P1_thresh{1,1,1} = 0; % Check this manually
% rsa_P1_thresh{3,1,2} = 0; % Check this manually

rsa_P1 = cellfun(@(x) (sum(x>r2t(0.05,nchoosek(ntrials,2)))./length(x))*100,rsa_P1_);
% rsa_P1 = cellfun(@(x1,x2) (sum(x1<x2.*ones(size(x1)))./length(x1))*100,rsa_P1_pvals,rsa_P1_thresh);
rsa_P1 = cellfun(@(x) (sum(x<rsa_P1_thresh)./length(x))*100,rsa_P1_pvals);
S_mat = squeeze(mean(rsa_P1));
S_err = squeeze(std(rsa_P1))./sqrt(3);
figure('Position',[0.5 0.5 400 250])
hold on
ngroups = size(S_mat, 1);
nbars = size(S_mat, 2);
bar(S_mat);
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
% legend({'Perceptual','Chemical','Mutual'})
% legend()
% Subject data points
c_s = {'r','g','b'}; % Data dots for subjects
for ii = 1:nanat % For bars for perceptual, chemical and combinations
    for jj = 1:3
        plot(x_m(:,ii),squeeze(rsa_P1(jj,ii,:)),c_s{jj})
    end
end
ylabel('Representational Similarity (%)')
% yline(r2t(0.05,sum(utl_mask2(:))));
yline(5);
savefig(fullfile(savepath,'ARC_RSA'))
print(fullfile(savepath,'ARC_RSA'),'-dpng')

%% Bar plots of Valence and Salience as a function of other dims
barplotter = false;

if barplotter
    chemer = false;
    valier = true;
    niter = 10000;
    nC = 20; % NUmber of chemical components

    figure('Position',[0.5 0.5 1280 360])
    hold on
    nS = 3;
    for ss = 1:nS
        subplot(1,3,ss)
        hold on
        Y = behavP.behav(ss).ratings(:,2);
        if ~valier;  Y = abs(Y); end
        if chemer
            X = behavC.behav(ss).ratings(:,1:nC); 
            xp = [1: nC]; % Labels
        else 
            X = behavP.behav(ss).ratings(:,[1 3:end]); 
            xp = behavP.behav(ss).percepts([1 3:end]);
        end
        figure()
        [beta_m, beta_err, p_values,~,res] = bootstrapRidgeres(X, Y, niter, 0.1);
        plotBarWithSignificance(beta_m, beta_err*1.96, p_values)
        % Note the p_value is approximate since it's not strictly LOOCV
        title(sprintf('S: %02d r: %.2f p: %.3f',ss,res,r2p(res,160)))
        if  ~chemer; xticks(1:17); xtickangle(90); xticklabels(xp); end
        ylabel('Ridge regression weights')
    end
    % savefig('Sal_perceptual2')
    % print('Sal_perceptual2','-dpng')
end



%% Normalize the images into MNI space
img_nrmr = false;
tic
for s = 1:3
    matlabbatch = [];
    dirs = {'C:\Data\NEMO\NEMO_01\imaging\nii\master_anat\meansNEMO_01_set1_sess1-0015-00001-000176-01.nii'
        'C:\Data\NEMO\NEMO_02\imaging\nii\anat\sNEMO02.nii'
        'C:\Data\NEMO\NEMO_04\imaging\nii\anat\sNEMO04.nii'};
    wdir = pwd;
%     f2 = dir(fullfile(wdir,sprintf('S%01d',s),'spmT_0006.nii'));
    f2 = dir(fullfile(wdir,sprintf('sfp_behav_s%02d',s),'rsa_vec_fless.nii'));

    files = {};
    for zz = 1:length(f2); files{zz,1} = fullfile(wdir,sprintf('sfp_behav_s%02d',s),f2(zz).name); end
    
    if img_nrmr
        matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = {dirs{s}};
        matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = files;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {'C:\Toolboxes\spm12\tpm\TPM.nii'};
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
        matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
            78 76 85];
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
        matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';
        spm_jobman('run', matlabbatch);
    end
end
toc
%% Image averager
img_avg = false;
ntrials = 4320;
thr = r2t(0.05,nchoosek(ntrials,2));
if img_avg
    outpurdir = 'C:\Data\SFP\neural';
    dirs = {'C:\Data\SFP\sfp_behav_s01'
            'C:\Data\SFP\sfp_behav_s02'
            'C:\Data\SFP\sfp_behav_s03'};
%           outpurdir = 'C:\Data\ARC\ARC\val_glm\val-sal\basic';
%     dirs = {'C:\Data\ARC\ARC\val_glm\val-sal\basic\S1'
%             'C:\Data\ARC\ARC\val_glm\val-sal\basic\S2'
%             'C:\Data\ARC\ARC\val_glm\val-sal\basic\S3'};
    files = {};
    for ss = 1:3
        f2 = dir(fullfile(dirs{ss},'wrsa_vec_fless.nii'));
        for zz = 1:length(f2); files{zz,ss} = fullfile(dirs{ss},f2(zz).name); end
    end
    
    for ff = 1:size(files,1)
        iminput = files(ff,:)';
        fname = files{ff,1}(end-8:end-4);
        matlabbatch = [];
        
        matlabbatch{1}.spm.util.imcalc.input = iminput;
        matlabbatch{1}.spm.util.imcalc.output = fname;
        matlabbatch{1}.spm.util.imcalc.outdir = {outpurdir};
%       matlabbatch{1}.spm.util.imcalc.expression = '(double(i1>0.00053)+double(i2>0.00053)+double(i3>0.00053))/3';
%         matlabbatch{1}.spm.util.imcalc.expression = '(-i1-i2-i3)/3';
         matlabbatch{1}.spm.util.imcalc.expression = '(i1+i2+i3)/3';
        matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
        matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
        matlabbatch{1}.spm.util.imcalc.options.mask = 0;
        matlabbatch{1}.spm.util.imcalc.options.interp = 1;
        matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
        spm_jobman('run', matlabbatch);
    end
end

% %% Behavioral stuff
% s = 3;
% C = behavC.behav(s).ratings(:,1:20);
% y = behav3(:,2);
% C3_scan = [C y];
% C3_scan(isnan(C3_scan))=0;
% writematrix(C3_scan)
% 
% % Misc rough work
% a = rand(160,1);
% for ii = 1:160; 
%     b = setdiff(1:160,ii); 
%     c(ii) = mean(a(b)); 
% end
% fastcorr(a,c')
