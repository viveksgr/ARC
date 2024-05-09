%% General Settings
root = 'C:\Work\ARC\ARC';
maskfile =  'ARC3_anatgw.nii';
fmaskfile = 'ARC3_fanatgw3_pos.nii';
fmasker = true;
binz = 7;
binzc = 0;
if mod(binz,2)==0; binzpart1 = binz/2; binzpart2 = binzpart1+1; else; binzpart1 = (binz+1)/2 ; binzpart2 = binzpart1; end

anat_names = {'PC','AMY','OFC','VMPFC'};
anat_masks = {'rwPC.nii','rwAmygdala.nii','rwofc.nii','rwvmpfc.nii'};
% anat_names = {'PC','AMY','Thal','frontal_sal','frontal_val'}; 
% anat_masks = {'rwPC.nii','rwAmygdala.nii','rwThal.nii','frontal_sal.nii','frontal_val.nii'};
nanat = length(anat_names);
medianize_behav = true;

rangenormer = false;
smoter = false;
pca_maker = false;
discretizer = false;
valsep = true;

% sess_l = cat(3,nchoosek([1 2 3],2),nchoosek([2 3 4],2),nchoosek([2 3
% 4],2),nchoosek([2 3 4],2)); % For sesswise
% load('C:\Data\NEMO\swampsunset.mat');
single_n = false; % Noisepool
single_c = true; % Cutoff from sign voxels
nodor = 160;

sess_l = repmat([0],1,2,3);
dirs = {fullfile(root,'\ARC01\mediation');
    fullfile(root,'\ARC01\mediation');
    fullfile(root,'\ARC01\mediation')};
behav = load(fullfile(root,'ARC','NEMO_perceptual2.mat'));
modelname = 'Decoding\ValSep_basic';
savepath = fullfile(root,modelname);
v_id = 2; % Index of vector for median splitting odors

% load(fullfile(statpath,'fir_cv.mat'))
fprintf('\n')

if ~valsep; rsa_P1 = zeros(3,nanat,3); else; rsa_P1 = zeros(3,nanat,4); end 
rsa_P1t = rsa_P1;
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
            modelmd2 = double(modelmd(~noisepool,:));
        else
            modelmd2 = double(modelmd);
        end
        modelmd2 = zscore(modelmd2);
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
                labels_val = discretize(behav_ratings_,binz);
                labels_sal = discretize(abs(behav_ratings_),binz);
                labels_combined =  (labels_val - 1) * binz + labels_sal;
            end
        else
            labels_val = behav_ratings_;
            labels_sal = abs(behav_ratings_);
        end

        if pca_maker
            [coeff,score,~,~,var] = pca(modelmd2');
            var = cumsum(var);
            nvox = min(size(score,2),100);
            modelmd2 = score(:,1:nvox )';
        end

        if smoter
            [neural_val,labels_val] = smote(modelmd2', [], 100, 'Class', labels_val);
            [neural_sal,labels_sal] = smote(modelmd2', [], 100, 'Class', labels_sal);
            [neural_val,labels_combined] = smote(modelmd2', [], 200, 'Class', labels_combined);
        else
            neural_val = modelmd2';
            neural_sal = modelmd2';
        end

    
        % 

        if ~valsep
            % % % % labels_val_perm = labels_val(randperm(length(labels_val)));
            [rsa_P1(s,ii,1),~,rsa_P1t(s,ii,1)] = Classify_Permute_VS2_regress(neural_val, labels_val, 4);

            neural_val_pos = neural_val(labels_val>binzc,:);
            neural_val_neg = neural_val(labels_val<binzc,:);

            labels_val_pos = labels_val(labels_val>binzc,:);
            labels_val_neg = labels_val(labels_val<binzc,:);

            [rsa_P1(s,ii,2),~,rsa_P1t(s,ii,2)] = Classify_Permute_VS2_regress( neural_val_pos, labels_val_pos, 4);
            [rsa_P1(s,ii,3),~,rsa_P1t(s,ii,3)] = Classify_Permute_VS2_regress( neural_val_neg, labels_val_neg, 4);
            
             % % Valsep
        else
            load('C:\Work\ARC\ARC\Decoding\t_score_mat.mat')
            assert(size(neural_val,2)==size(t_score_mat{s,ii},1))
            thr = tinv(0.975,size(t_score_mat{s,ii},1));
            pospop = and(t_score_mat{s,ii}(:,1)>thr,t_score_mat{s,ii}(:,2)<thr);
            negpop = and(t_score_mat{s,ii}(:,2)>thr,t_score_mat{s,ii}(:,1)<thr);
            neural_val_pos = neural_val(labels_val>binzc,pospop); % Pos pop coding pos
            neural_val_neg = neural_val(labels_val<binzc,negpop); % Neg pop coding neg
            neural_val_pos_cr = neural_val(labels_val>binzc,negpop); % Neg pop coding pos
            neural_val_neg_cr = neural_val(labels_val<binzc,pospop); % Pos pop coding neg

            labels_val_pos = labels_val(labels_val>binzc);
            labels_val_neg = labels_val(labels_val<binzc);

            [rsa_P1(s,ii,1),~,rsa_P1t(s,ii,1)] = Classify_Permute_VS2_regress(neural_val_pos, labels_val_pos, 10);
            [rsa_P1(s,ii,2),~,rsa_P1t(s,ii,2)] = Classify_Permute_VS2_regress(neural_val_neg, labels_val_neg, 10);
            [rsa_P1(s,ii,3),~,rsa_P1t(s,ii,3)] = Classify_Permute_VS2_regress(neural_val_pos_cr, labels_val_pos, 10);
            [rsa_P1(s,ii,4),~,rsa_P1t(s,ii,4)] = Classify_Permute_VS2_regress(neural_val_neg_cr, labels_val_neg, 10);
            
            

            % mdl = svmtrain( trainY,  trainX, '-s 3 -t 2 -c 1 -q');
            % % eval( ['mdl = svmtrain( l1( :, 1), mat( l1( :, 2), :), '' ', opt.svm, ''');']);
            % predictions = svmpredict(  testY , testX ,  mdl, ' -q ');
            % rsa_P1(s,ii,3) = corr(predictions,testY);
            % 
            % mdl = svmtrain( testY,  testX, '-s 3 -t 2 -c 1 -q');
            % predictions= svmpredict(  trainY , trainX ,  mdl, ' -q ');
            % rsa_P1(s,ii,4) = corr(predictions,trainY);
        end


        % % Combined classification
        % [~,pred] = Classify_Permute_VS2(modelmd2', labels_combined, 4);
        % [pred_val, pred_sal] = ARC_extractLabels(pred,binz);
        % [act_val,act_sal] = ARC_extractLabels(labels_combined,binz);
        % 
        % rsa_P1(s,ii,1) = sum(pred_val==act_val)/length(act_val);
        % rsa_P1(s,ii,2) = sum(pred_sal==act_sal)/length(act_sal);

        % % % Sanity checks
        % modelmd_binned = ARC_binAndTransform(modelmd2, behav_ratings_, binz, [ms1 ms2]);
        % modelmd_corrcoef = corrcoef(modelmd_binned);
        % val_sc = linspace(-1,1,binz);
        % sal_sc = abs(val_sc);
        % val_mat = 1-abs(val_sc-val_sc');
        % sal_mat = 1-abs(sal_sc-sal_sc');
        % utl_mask2 = logical(triu(ones(length(val_mat)),1)); % All possible odors
        % % 
        % % fastcorr(sal_mat(utl_mask2),modelmd_corrcoef(utl_mask2))
        % % fastcorr(val_mat(utl_mask2),modelmd_corrcoef(utl_mask2))
        % % fastcorr(val_mat(utl_mask2),sal_mat(utl_mask2))
        % 
        % [w,t] = ARC_multicomputeWeights_tsc( [val_mat(utl_mask2) sal_mat(utl_mask2)],modelmd_corrcoef(utl_mask2));
    end
end

plotsw = false;
if plotsw; rsa_plt = rsa_P1t; else; rsa_plt = rsa_P1; end

mkdir(savepath)
S_mat = squeeze(mean(rsa_plt));
S_err = squeeze(std(rsa_plt))./sqrt(3);
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
        plot(x_m(:,ii),squeeze(rsa_plt(jj,ii,:)),c_s{jj})
    end
end
if  plotsw
    ylabel('Decoding accuracy (t)')
else
    ylabel('Decoding accuracy (r)')
end
if ~valsep; legend({'Val','Val+','Val-'}); else; legend({'Val+ from N+','Val- from N-','Val+ from N-','Val- from N+'}); end
if  plotsw
yline(1.65,'HandleVisibility','off')
yline(-1.65,'HandleVisibility','off') 
end
% yline(r2t(0.025,4320),'HandleVisibility','off')
% yline(-r2t(0.025,4320),'HandleVisibility','off')
% legend({'Valence','Salience'})
clear modelmd_ modelmd modelmd2 S1_omat_vals S2_omat_vals unity M_anat M_sal M_val 
p_values = arrayfun(@(r_eff) ARC_computePValueOneTailed(r_eff, binz, length(behav_ratings_)), rsa_P1);
savefig(fullfile(savepath,'ARC_decoding'))

save(fullfile(savepath,'ARC_decoding'))