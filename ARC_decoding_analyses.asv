%% General Settings

% Check fmasker
tic
root = 'D:\Work\ARC\ARC'; %'D:\Work\ARC\Scripts\'; %
maskfile =  'ARC3_anatgw.nii';
fmaskfile = 'ARC3_fanatgw3_pos.nii';
fmasker = false; % Don't use functional masking
binz = 7;
binzc = 4;
if mod(binz,2)==0; binzpart1 = binz/2; binzpart2 = binzpart1+1; else; binzpart1 = (binz+1)/2 ; binzpart2 = binzpart1; end
 
% anat_names = {'Insula','Hipp','DLPFC','A1','wm'};
% anat_masks = {'rwinsula.nii','rwHipp.nii','rwDLPFC.nii','rwAud.nii','rwm_main.nii'};

anat_names = {'PC','AMY','OFC','VMPFC'};
anat_masks = {'rwPC.nii','rwAmygdala.nii','rwofc.nii','rwvmpfc.nii'};

% anat_names = {'11m','14m','47m','47o','ar11','ar13','FPI','FPM'};
% anat_masks = {'rwar_03_24.nii','rwar_04_25.nii','rwar_08_29.nii','rwar_09_30.nii','rwar_10_31.nii','rwar_11_32.nii','rwar_15_36.nii','rwar_16_37.nii'};
% anat_names = {'Clust1','Clust2','Clust3','Clust4','Clust5','Clust6'};
% anat_masks = {'rwk_clus1.nii','rwk_clus2.nii','rwk_clus3.nii','rwk_clus4.nii','rwk_clus5.nii','rwk_clus6.nii'};

ar_exOFC = false;
nanat = length(anat_names);
incneut = false;
medianize_behav = true;
shuffler = false;
rangenormer = false;

pca_maker = false;
discretizer = false;
sesswise = false;
switcher = 'Domainpart'; % 'Basic'; 'Domainpart'; 'Neuralpart_mutual'; 'Crossdec';'CrossdecNeural'
savename = 'final_f_decoding/Basic_rep';

single_n = false; % Noisepool
single_c = false; % Cutoff from sign voxels
z_scorer = false;
nodor = 160;

sz_ctrl = false;
intens_reg = false;
wt_adj = false;

sess_l = repmat([0],1,2,3);


dirs = {fullfile(root,'\ARC01\mediation');
    fullfile(root,'\ARC02\mediation');
    fullfile(root,'\ARC03\mediation')};
behav = load(fullfile(root,'ARC','NEMO_perceptual2.mat'));
modelname = fullfile('Decoding',switcher);
savepath = fullfile(root,savename);
v_id = 2; % Index of vector for median splitting odors

% load(fullfile(statpath,'fir_cv.mat'))
fprintf('\n')

rsa_P1 = zeros(3,nanat,2);
rsa_P1t = rsa_P1;
dfs = rsa_P1;
rs = zeros(3,1);
hold on
% Subject - index
kk = 1;

if fmasker
    % ffile ='C:\Work\ARC\ARC\RSA\controls\main\ARC_RSA.mat';
    % ffile = 'C:\Work\ARC\ARC\RSA\Basic\ARC_RSA.mat';
    ffile = 'D:\Work\ARC\ARC\RSA\valsep_corrected2\ARC_RSA.mat';
else
    % ffile = 'C:\Work\ARC\ARC\RSA\unmasked\valsep\ARC_RSA.mat';
    % ffile = 'C:\Work\ARC\ARC\RSA\unmasked\basic_med\ARC_RSA.mat';
    ffile = 'D:\Work\ARC\ARC\results\basic_main_100\ARC_RSA.mat';
end


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
    m_3d = {};
    for ii = 1:length(anat_masks)
        m1 = spm_read_vols(spm_vol(fullfile(anatdir,anat_masks{ii})));
        m1(isnan(m1))=0;
        m1(m1<=0.01)=0;
        m1(m1>0) = 1;
        m_3d{ii} = and(m1,mask);
        m1 = m1(mask);
        masks_set(:,ii)=m1(fmask_1d);
        fprintf('area count: %04d\n',sum(m1(fmask_1d)))
        masks_set_cell{ii} = fmask_1d(logical(m1));
      
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

    % Shufle sanity check
    if shuffler
        behav_ratings_ = behav_ratings_(randperm(length(behav_ratings_)));
    end
    sum_vox = [];
    for ii = 1:length(anat_names)
        fprintf('area:%02d\n',ii)
        % if and(s==2,ii==7)
        %     'beep'
        % end
        if ar_exOFC
            load(fullfile(anatdir,'ar_exmats','ar_ex_mats.mat'));
            % --- select rows whose voxels are inside m3, preserving out.idx order
            uIdx = out.idx(:);             
            % union indices that define row order of out.X
            m3Idx = find(m_3d{ii});
            keep = ismember(uIdx, m3Idx);    % logical over rows of out.X
             modelmd = out.X(keep, :);
            sum_vox(ii) = sum(keep);
        else
        % modelmd_ = load(fullfile(anatdir,'sesswise',anat_names{ii},'TR_01','TYPED_FITHRF_GLMDENOISE_RR.mat'),'modelmd','noisepool');
         modelmd_ = load(fullfile(anatdir,anat_names{ii},'TYPED_FITHRF_GLMDENOISE_RR.mat'),'modelmd','noisepool');
      
        modelmd = squeeze(modelmd_.modelmd);
        noisepool = modelmd_.noisepool;
        end

        if single_c
            modelmd = modelmd(masks_set_cell{ii},:);
            noisepool = noisepool(masks_set_cell{ii});
        end

        if single_n
            modelmd2 = double(modelmd(~noisepool,:));
        else
            modelmd2 = double(modelmd);
        end
        if z_scorer
            modelmd2 = zscore(modelmd2);
        end

        fprintf('size:%02d\n',size(modelmd,1))

        [r1,~] = find(isnan(modelmd2));
        modelmd2(r1,:) = [];

        if discretizer
            if rangenormer
                edges_val = quantile(behav_ratings_,linspace(0,1,binz+1));
                edges_sal = quantile(abs(behav_ratings_),linspace(0,1,binz+1));
                labels_val = discretize(behav_ratings_,edges_val);
                labels_sal = discretize(abs(behav_ratings_),edges_sal);
            else
                labels_int = behav.behav(s).ratings(group_vec,1);
                if intens_reg
                    labels_val = discretize(regressmeout(behav_ratings_',labels_int')',binz);
                    labels_sal = discretize(regressmeout(abs(behav_ratings_)',labels_int')',binz);
                else
                    labels_val = discretize(behav_ratings_,binz);
                    labels_sal = discretize(abs(behav_ratings_),binz);
                end
                labels_combined =  (labels_val - 1) * binz + labels_sal;
            end
        else
            
            if intens_reg
                  labels_int = behav.behav(s).ratings(group_vec,1);
                labels_val = regressmeout(behav_ratings_',labels_int')';
                labels_sal = regressmeout(abs(behav_ratings_)',labels_int')';
            else
                labels_val = behav_ratings_;
                labels_sal = abs(behav_ratings_);
            end
        end

        if pca_maker
            [coeff,score,~,~,var] = pca(modelmd2');
            var = cumsum(var);
            nvox = min(size(score,2),100);
            modelmd2 = score(:,1:nvox )';
        end

        if sesswise
            R = ARC_extract_roi_trials_travg_dec(anatdir,anat_masks,ii,s);
            neural_val = R.betatrials';
            behav_ratings_sess = behav_ratings(R.group_vec);
            labels_val = discretize(behav_ratings_sess,binz);
            labels_sal= discretize(abs(behav_ratings_sess),binz);
        else
            neural_val = modelmd2';
        end

        switch switcher
            case 'Basic'
                if sz_ctrl
                    nperm = 100;
                    tempmat1 = zeros(nperm,1);    tempmat2 = zeros(nperm,1);
                    for pp = 1:nperm
                        vox_id =datasample(1:size(neural_val,2),nperm);
                        tempmat1(pp) = ARC_regress_normed(neural_val(:,vox_id), labels_val, 2,1);
                        tempmat2(pp) = ARC_regress_normed(neural_val(:,vox_id), labels_sal, 2,1);
                    end
                    rsa_P1(s,ii,1) = mean( tempmat1);
                    rsa_P1(s,ii,2) = mean( tempmat2);
                    [~,rsa_P1t(s,ii,1)] = ARC_r2t(mean( tempmat1),length(labels_val));
                    [~,rsa_P1t(s,ii,2)] = ARC_r2t( mean( tempmat2),length(labels_sal));
                else
                    if wt_adj
                        nperm = 20;
                        tempmat1 = zeros(nperm,1);    tempmat2 = zeros(nperm,1);
                        for pp = 1:nperm

                        [~,idx_val ] = datasample(labels_val,1000);
                        [~,idx_sal] =datasample(labels_sal,1000);
                           tempmat1(pp) = ARC_regress_normed(neural_val(idx_val,:), labels_val(idx_val), 2,1);
                            tempmat2(pp) = ARC_regress_normed(neural_val(idx_sal,:), labels_sal(idx_sal), 2,1);
                        end
                        rsa_P1(s,ii,1) = mean( tempmat1);
                        rsa_P1(s,ii,2) = mean( tempmat2);
                        [~,rsa_P1t(s,ii,1)] = ARC_r2t(mean( tempmat1),length(labels_val));
                        [~,rsa_P1t(s,ii,2)] = ARC_r2t( mean( tempmat2),length(labels_sal));
                    else
                        [rsa_P1(s,ii,1),~,rsa_P1t(s,ii,1)] = ARC_regress_normed(neural_val, labels_val, 2,1);
                        [rsa_P1(s,ii,2),~,rsa_P1t(s,ii,2)] = ARC_regress_normed(neural_val, labels_sal, 2,1);
                    end
                end

                legender = {'Val','Sal'};
                dfs(s,ii,1)=length(labels_val);
                dfs(s,ii,2)=length(labels_sal);
            case 'Domainpart'
                % % % % labels_val_perm = labels_val(randperm(length(labels_val)));
                if ~discretizer; binzc = 0; end
                neural_val_pos = neural_val(labels_val>=binzc,:);
                neural_val_neg = neural_val(labels_val<=binzc,:);
                labels_val_pos = labels_val(labels_val>=binzc,:); 
                labels_val_neg = labels_val(labels_val<=binzc,:);

                if sz_ctrl
                    nperm = 100;
                    tempmat1 = zeros(nperm,1);    tempmat2 = zeros(nperm,1);
                    for pp = 1:nperm
                        vox_id =datasample(1:size(neural_val,2),100);
                        tempmat1(pp) = ARC_regress_normed( neural_val_pos(:,vox_id), labels_val_pos, 2,1);
                        tempmat2(pp) = ARC_regress_normed( neural_val_neg(:,vox_id), labels_val_neg, 2,1);
                    end
                    rsa_P1(s,ii,1) = mean( tempmat1);
                    rsa_P1(s,ii,2) = mean( tempmat2);
                    [~,rsa_P1t(s,ii,1)] = ARC_r2t(mean( tempmat1),length(labels_val_pos));
                    [~,rsa_P1t(s,ii,2)] = ARC_r2t( mean( tempmat2),length(labels_val_neg));

                    dfs(s,ii,1)=length(labels_val_neg);
                    dfs(s,ii,2)=length(labels_val_pos);
                else
                    if wt_adj

                        [~,idx_valp] = datasample(labels_val_pos,1000);
                        [~,idx_valn] =datasample(labels_val_neg,1000);

                        [rsa_P1(s,ii,1),~,rsa_P1t(s,ii,1)] = ARC_regress_normed(neural_val_pos(idx_valp,:), labels_val_pos(idx_valp), 2,1);
                        [rsa_P1(s,ii,2),~,rsa_P1t(s,ii,2)] = ARC_regress_normed(neural_val_neg(idx_valn,:), labels_val_neg(idx_valn), 2,1);
                   

                        dfs(s,ii,1)=length(idx_valp);
                        dfs(s,ii,2)=length(idx_valn);


                    else

                        [rsa_P1(s,ii,1),~,rsa_P1t(s,ii,1)] = ARC_regress_normed(neural_val_pos, labels_val_pos, 2,1);
                        [rsa_P1(s,ii,2),~,rsa_P1t(s,ii,2)] = ARC_regress_normed(neural_val_neg, labels_val_neg, 2,1);
                        
                        dfs(s,ii,1)=length(labels_val_neg);
                        dfs(s,ii,2)=length(labels_val_pos);
                    end
                end
                legender = {'Val+','Val-'};

            case 'Crossdec'
                if discretizer; bcn = binzc; else; bcn = 0; end

                if incneut
                    labneg = selectLess4AndHalf4(labels_val,bcn);
                    neural_val_pos = neural_val(~labneg ,:);
                    neural_val_neg = neural_val(labneg ,:);

                    labels_val_pos = labels_val(~labneg ,:);
                    labels_val_neg = labels_val(labneg ,:);

                else
                    neural_val_pos = neural_val(labels_val>binzc,:);
                    neural_val_neg = neural_val(labels_val<binzc,:);

                    labels_val_pos = labels_val(labels_val>binzc,:);
                    labels_val_neg = labels_val(labels_val<binzc,:);
                end


                [ rsa_P1(s,ii,1),~,rsa_P1t(s,ii,1)] = ARC_regress_wrapper(neural_val_pos,labels_val_pos, neural_val_neg, labels_val_neg, 10,1,false);
                [ rsa_P1(s,ii,2),~,rsa_P1t(s,ii,2)] = ARC_regress_wrapper(neural_val_neg, labels_val_neg,neural_val_pos,labels_val_pos, 10,1,false);
         
                legender = {'Val- from Val+','Val- from Val+'};

                dfs(s,ii,1)=length(labels_val_neg);
                dfs(s,ii,2)=length(labels_val_pos);

          
            otherwise
                error('Decoding type unspecified...')
        end
    end
end


mkdir(savepath)
% savepath = pwd;
rsa_copy = rsa_P1;
p_values_3dt = ARC_decoding_pvals(rsa_P1, dfs);
% rsa_P1(vpop<10) = nan;
ARC_barplot_sig(rsa_P1,p_values_3dt,true,false)
xticks(1:nanat)
xticklabels(anat_names);
% ylim([-0.1 0.4])
ylabel('Prediction (r)')
% legend({'Val+ -> Val-','Val- -> Val+'})
legend(legender,'Location','southwest')
% ylim([-0.1 0.4])
savefig(fullfile(savepath,'ARC_decoding'))
print(gcf,'-vector','-dsvg',[fullfile(savepath,'ARC_decoding_diff'),'.svg']) % svg
print([fullfile(savepath,'ARC_decoding_diff')],'-dpng') % svg
% print(gcf,'-vector','-dsvg',[fullfile(savepath,'ARC_decoding'),'.svg']) % svg
% p_values_3d2 = ARC_combinePValues3D(rsa_P1t, rsa_P1, dfs)
% p_values_3dz = ARC_decoding_pvals_zsc(rsa_P1, dfs)
% savepath = pwd;-
save(fullfile(savepath,'ARC_decoding'))
toc

if strcmp(switcher,'Crossdec')
    beta_mat  = tanh(mean(atanh(rsa_P1),3));
    % p_vals = max(p_values_3dt ,[],2);
    p_values_3dt = ARC_decoding_pvals(  beta_mat, mean(dfs,3));
    plot_roi_bars_with_subjects(beta_mat, p_vals, ...
    'RoiLabels',anat_names, ...
    'YLabel','Accuracy (r)', ...
    'Title','Cross Decoding', ...
    'ErrType','sem', ...
    'StarColor',[0.8 0 0], 'StarFontSize',12, 'ErrMult',1.6);
end

% Single r plot
if strcmp(switcher,'Domainpart')
    % [beta_mat, p_vals] = ARC_diff2conds_dfweighted(rsa_P1, dfs);
    [beta_mat, p_vals] = ARC_diff2conds_lessConservative(rsa_P1, dfs, ...
                              'Rho',0.3, 'Method','fixed');
    plot_roi_bars_with_subjects(beta_mat, p_vals, ...
    'RoiLabels',anat_names, ...
    'YLabel','Accuracy (r)', ...
    'Title','beta_va', ...
    'ErrType','sem', ...
    'StarColor',[0.8 0 0], 'StarFontSize',12, 'ErrMult',1.6);
    % ylim([-4 4])
    print(gcf,'-vector','-dsvg',[fullfile(savepath,'feldspar'),'.svg']) % svg

end

