%% General Settings

% Check fmasker
tic
root = 'C:\Work\ARC\ARC';
maskfile =  'ARC3_anatgw.nii';
fmaskfile = 'ARC3_fanatgw3_pos.nii';
fmasker = false;
binz = 7;
binzc = 4;
if mod(binz,2)==0; binzpart1 = binz/2; binzpart2 = binzpart1+1; else; binzpart1 = (binz+1)/2 ; binzpart2 = binzpart1; end

% anat_names = {'Insula','Hipp','DLPFC','A1','wm'};
% anat_masks = {'rwinsula.nii','rwHipp.nii','rwDLPFC.nii','rwAud.nii','rwm_main.nii'};

anat_names = {'PC','AMY','OFC','VMPFC'};
anat_masks = {'rwPC.nii','rwAmygdala.nii','rwofc.nii','rwvmpfc.nii'};

nanat = length(anat_names);
medianize_behav = true;
shuffler = false;
rangenormer = false;
smoter = false;
pca_maker = false;
discretizer = true;
switcher = 'Domainpart'; % 'Basic'; 'Domainpart'; 'Neuralpart_mutual'; 'Crossdec';'CrossdecNeural'
popname = 'pos'; % 'pos' 'neg' or 'mut'
savename = 'RSA_correct/Domainpart_main_woneut';
% sess_l = cat(3,nchoosek([1 2 3],2),nchoosek([2 3 4],2),nchoosek([2 3
% 4],2),nchoosek([2 3 4],2)); % For sesswise
% load('C:\Data\NEMO\swampsunset.mat');
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
savepath = fullfile(root,'Decoding',savename);
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
    ffile = 'C:\Work\ARC\ARC\RSA\valsep_corrected2\ARC_RSA.mat';
else
    % ffile = 'C:\Work\ARC\ARC\RSA\unmasked\valsep\ARC_RSA.mat';
    % ffile = 'C:\Work\ARC\ARC\RSA\unmasked\basic_med\ARC_RSA.mat';
    ffile = 'C:\Work\ARC\ARC\results\basic_main_100\ARC_RSA.mat';
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


    % Shufle sanity check
    if shuffler
        behav_ratings_ = behav_ratings_(randperm(length(behav_ratings_)));
    end
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
            labels_val = behav_ratings_;
            labels_sal = abs(behav_ratings_);
            labels_int = behav.behav(s).ratings(group_vec,1);
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

        switch switcher
            case 'Basic'
                if sz_ctrl
                    nperm = 100;
                    tempmat1 = zeros(nperm,1);    tempmat2 = zeros(nperm,1);
                    for pp = 1:nperm
                        vox_id =datasample(1:size(neural_val,2),nperm);
                        tempmat1(pp) = ARC_regress_nested2_normed(neural_val(:,vox_id), labels_val, 2,1);
                        tempmat2(pp) = ARC_regress_nested2_normed(neural_val(:,vox_id), labels_sal, 2,1);
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
                            idx_val = ARC_balanceClasses(labels_val, binz);
                            idx_sal = ARC_balanceClasses(labels_sal, binz);
                            tempmat1(pp) = ARC_regress_nested2_normed(neural_val(idx_val,:), labels_val(idx_val), 2,1);
                            tempmat2(pp) = ARC_regress_nested2_normed(neural_val(idx_sal,:), labels_sal(idx_sal), 2,1);
                        end
                        rsa_P1(s,ii,1) = mean( tempmat1);
                        rsa_P1(s,ii,2) = mean( tempmat2);
                        [~,rsa_P1t(s,ii,1)] = ARC_r2t(mean( tempmat1),length(labels_val));
                        [~,rsa_P1t(s,ii,2)] = ARC_r2t( mean( tempmat2),length(labels_sal));
                    else
                        % [rsa_P1(s,ii,1),~,rsa_P1t(s,ii,1)] = ARC_regress_nested2_normed(neural_val, labels_val, 2,1);
                        % [rsa_P1(s,ii,2),~,rsa_P1t(s,ii,2)] = ARC_regress_nested2_normed(neural_val, labels_sal, 2,1);
                        % if and(ii==1, s==1)
                        % 'test'
                        % end
                         [rsa_P1(s,ii,1),~,rsa_P1t(s,ii,1)] = ARC_regress_normed(neural_val, labels_val, 2,1);
                        [rsa_P1(s,ii,2),~,rsa_P1t(s,ii,2)] = ARC_regress_normed(neural_val, labels_sal, 2,1);

                        % [rsa_P1(s,ii,1),~,rsa_P1t(s,ii,1)] = Classify_Permute_VS2_regress(neural_val, labels_val, 4);
                        % [rsa_P1(s,ii,2),~,rsa_P1t(s,ii,2)] = Classify_Permute_VS2_regress(neural_val, labels_sal, 4);
                    end
                end

                % [pred_val, pred_sal] = ARC_extractLabels(pred,binz);
                % [act_val,act_sal] = ARC_extractLabels(labels_combined,binz);
                % %
                % rsa_P1(s,ii,1) = sum(pred_val==act_val)/length(act_val);
                % rsa_P1(s,ii,2) = sum(pred_sal==act_sal)/length(act_sal);
                %
                % p_accu = fastcorr(predictions,grp);
                % [~,p_accut] = ARC_r2t(p_accu,length(grp));

                legender = {'Val','Sal'};
                dfs(s,ii,1)=length(labels_val);
                dfs(s,ii,2)=length(labels_sal);
            case 'Domainpart'
                % % % % labels_val_perm = labels_val(randperm(length(labels_val)));
                neural_val_pos = neural_val(labels_val>binzc,:);
                neural_val_neg = neural_val(labels_val<binzc,:);

                labels_val_pos = labels_val(labels_val>binzc,:);
                labels_val_neg = labels_val(labels_val<binzc,:);

                if sz_ctrl
                    nperm = 100;
                    tempmat1 = zeros(nperm,1);    tempmat2 = zeros(nperm,1);
                    for pp = 1:nperm
                        vox_id =datasample(1:size(neural_val,2),100);
                        tempmat1(pp) = ARC_regress_nested2_normed( neural_val_pos(:,vox_id), labels_val_pos, 2,1);
                        tempmat2(pp) = ARC_regress_nested2_normed( neural_val_neg(:,vox_id), labels_val_neg, 2,1);
                    end
                    rsa_P1(s,ii,1) = mean( tempmat1);
                    rsa_P1(s,ii,2) = mean( tempmat2);
                    [~,rsa_P1t(s,ii,1)] = ARC_r2t(mean( tempmat1),length(labels_val_pos));
                    [~,rsa_P1t(s,ii,2)] = ARC_r2t( mean( tempmat2),length(labels_val_neg));

                    dfs(s,ii,1)=length(labels_val_neg);
                    dfs(s,ii,2)=length(labels_val_pos);
                else
                    if wt_adj
                        labels_val_pos = labels_val_pos-3;
                        nperm = 20;
                        tempmat1 = zeros(nperm,1);    tempmat2 = zeros(nperm,1);
                        for pp = 1:nperm
                            idx_valp = ARC_balanceClasses(labels_val_pos, length(unique(labels_val_pos)),200);
                            idx_valn = ARC_balanceClasses(labels_val_neg, length(unique(labels_val_neg)),200);
                            tempmat1(pp) = ARC_regress_nested2_normed(neural_val_pos(idx_valp,:), labels_val_pos(idx_valp), 2,1);
                            tempmat2(pp) = ARC_regress_nested2_normed(neural_val_neg(idx_valn,:), labels_val_neg(idx_valn), 2,1);
                        end
                        rsa_P1(s,ii,1) = mean( tempmat1);
                        rsa_P1(s,ii,2) = mean( tempmat2);
                        [~,rsa_P1t(s,ii,1)] = ARC_r2t(mean( tempmat1),length(labels_val_pos));
                        [~,rsa_P1t(s,ii,2)] = ARC_r2t( mean( tempmat2),length(labels_val_neg));

                        dfs(s,ii,1)=length(idx_valp);
                        dfs(s,ii,2)=length(idx_valn);
                    else
                        % [rsa_P1(s,ii,1),~,rsa_P1t(s,ii,1)] = Classify_Permute_VS2_regress( neural_val_pos, labels_val_pos, 4);
                        % [rsa_P1(s,ii,2),~,rsa_P1t(s,ii,2)] = Classify_Permute_VS2_regress( neural_val_neg, labels_val_neg, 4);
                        % 'test'
                        %  if and(ii==1, s==1)
                        % 'test'
                        %  end

                        [rsa_P1(s,ii,1),~,rsa_P1t(s,ii,1)] = ARC_regress_normed(neural_val_pos, labels_val_pos, 4,1);
                        [rsa_P1(s,ii,2),~,rsa_P1t(s,ii,2)] = ARC_regress_normed(neural_val_neg, labels_val_neg, 4,1);
                        % 
                        % [rsa_P1(s,ii,1),~,rsa_P1t(s,ii,1)] = ARC_regress_nested2_normed(neural_val_pos, labels_val_pos, 4,1);
                        % [rsa_P1(s,ii,2),~,rsa_P1t(s,ii,2)] = ARC_regress_nested2_normed(neural_val_neg, labels_val_neg, 4,1);

                        dfs(s,ii,1)=length(labels_val_neg);
                        dfs(s,ii,2)=length(labels_val_pos);
                    end
                end
                legender = {'Val+','Val-'};

                % % Valsep
            case 'Neuralpart_mutual'
                if and(s==3,ii==3)
                    map = 'beep';
                end
                load(ffile,'t_score_mat','thr_fdranat','w_score_mat')
                assert(size(neural_val,2)==size(t_score_mat{s,ii},1))

                % thr = tinv(0.9,size(t_score_mat{s,ii},1));
                % thr1 = thr_fdranat(s,1,ii);
                % thr2 = thr_fdranat(s,2,ii);
                % thr1 = 1.98;
                % thr2 = 1.98;
                % % pospop = and(t_score_mat{s,ii}(:,1)>thr1,t_score_mat{s,ii}(:,2)<thr);
                % % negpop = and(t_score_mat{s,ii}(:,2)>thr2,t_score_mat{s,ii}(:,1)<thr);
                % pospop = and(t_score_mat{s,ii}(:,1)>thr1,t_score_mat{s,ii}(:,2)>thr2);
                % negpop = and(t_score_mat{s,ii}(:,2)>thr1,t_score_mat{s,ii}(:,1)>thr2);

                thr = 0.075;
                thr1 = 0.075;
                thr2 = 0.075;
                pospop = and(w_score_mat{s,ii}(:,1)>thr1,w_score_mat{s,ii}(:,2)<thr);
                negpop = and(w_score_mat{s,ii}(:,2)>thr2,w_score_mat{s,ii}(:,1)<thr);
                mutpop = and(w_score_mat{s,ii}(:,1)>thr1,w_score_mat{s,ii}(:,2)>thr2);

                vpop(s,ii,1) = sum(pospop);
                vpop(s,ii,2) = sum(negpop);

                neural_val_pos = neural_val(labels_val>binzc,pospop); % Pos pop coding pos
                neural_val_neg = neural_val(labels_val<binzc,negpop); % Neg pop coding neg
                neural_val_pos_cr = neural_val(labels_val>binzc,negpop); % Neg pop coding pos
                neural_val_neg_cr = neural_val(labels_val<binzc,pospop); % Pos pop coding neg
                neural_val_mutpos = neural_val(labels_val>binzc,mutpop); % Pos pop coding pos
                neural_val_mutneg = neural_val(labels_val<binzc,mutpop); % Neg pop coding neg


                labels_val_pos = labels_val(labels_val>binzc);
                labels_val_neg = labels_val(labels_val<binzc);

                switch popname
                    case 'mut'
                        % Mutpop
                        [rsa_P1(s,ii,1),~,rsa_P1t(s,ii,1)] = ARC_regress_nested2_normed(neural_val_mutpos, labels_val_pos, 2,1);
                        [rsa_P1(s,ii,2),~,rsa_P1t(s,ii,2)] = ARC_regress_nested2_normed(neural_val_mutneg, labels_val_neg, 2,1);
                         legender = {'Val+ from N0','Val- from N0'};
                    case 'pos'

                        % Pospop
                        [rsa_P1(s,ii,1),~,rsa_P1t(s,ii,1)] = ARC_regress_nested2_normed(neural_val_pos, labels_val_pos, 2,1);
                        [rsa_P1(s,ii,2),~,rsa_P1t(s,ii,2)] = ARC_regress_nested2_normed(neural_val_neg_cr, labels_val_neg, 2,1);
                         legender = {'Val+ from N+','Val- from N+'};
                    case 'neg'
                        % Negpop
                        [rsa_P1(s,ii,1),~,rsa_P1t(s,ii,1)] = ARC_regress_nested2_normed(neural_val_pos_cr, labels_val_pos, 2,1);
                        [rsa_P1(s,ii,2),~,rsa_P1t(s,ii,2)] = ARC_regress_nested2_normed(neural_val_neg, labels_val_neg, 2,1);
                         legender = {'Val+ from N-','Val- from N-'};
                    otherwise
                        error('Wrong pop chosen')
                end              
                dfs(s,ii,1)=length(labels_val_pos);
                dfs(s,ii,2)=length(labels_val_neg);

            case 'Crossdec'
                labneg = selectLess4AndHalf4(labels_val);



                % neural_val_pos = neural_val(labels_val>binzc,:);
                % neural_val_neg = neural_val(labels_val<binzc,:);
                % 
                % labels_val_pos = labels_val(labels_val>binzc,:);
                % labels_val_neg = labels_val(labels_val<binzc,:);


                neural_val_pos = neural_val(~labneg ,:);
                neural_val_neg = neural_val(labneg ,:);

                labels_val_pos = labels_val(~labneg ,:);
                labels_val_neg = labels_val(labneg ,:);

                [ rsa_P1(s,ii,1),~,rsa_P1t(s,ii,1)] = ARC_regress_wrapper(neural_val_pos,labels_val_pos, neural_val_neg, labels_val_neg, 4,1,false);
                [ rsa_P1(s,ii,2),~,rsa_P1t(s,ii,2)] = ARC_regress_wrapper(neural_val_neg, labels_val_neg,neural_val_pos,labels_val_pos, 4,1,false);

                % % mdl = svmtrain(labels_val_pos, neural_val_pos,   '-s 3 -t 2 -c 1 -q');
                % mdl = ARC_regress_nested_mdl(neural_val_pos, labels_val_pos, 2);
                % predictions = svmpredict(  labels_val_neg , neural_val_neg ,  mdl, ' -q ');
                %
                % rsa_P1(s,ii,1) = fastcorr(predictions,labels_val_neg);
                % [~,rsa_P1t(s,ii,1)] = ARC_r2t(rsa_P1(s,ii,1),length(labels_val_neg));
                %
                % % mdl = svmtrain(labels_val_neg, neural_val_neg,   '-s 3 -t 2 -c 1 -q');
                % mdl = ARC_regress_nested_mdl(neural_val_neg, labels_val_neg, 2);
                % predictions = svmpredict(  labels_val_pos , neural_val_pos ,  mdl, ' -q ');
                %
                % rsa_P1(s,ii,2) = fastcorr(predictions,labels_val_pos);
                % [~,rsa_P1t(s,ii,2)] = ARC_r2t(rsa_P1(s,ii,2),length(labels_val_pos));

                legender = {'Val- from Val+','Val- from Val+'};

                dfs(s,ii,1)=length(labels_val_neg);
                dfs(s,ii,2)=length(labels_val_pos);

            case 'CrossdecNeural'
                load(ffile,'t_score_mat','thr_fdranat','w_score_mat')
                assert(size(neural_val,2)==size(t_score_mat{s,ii},1))
                thr = tinv(0.90,size(t_score_mat{s,ii},1));

                % thr1 = thr_fdranat(s,1,ii);
                % thr2 = thr_fdranat(s,2,ii);

                % thr1 = 1.65;
                % thr2 = 1.65;
                % thr = 1.65;
                % pospop = and(t_score_mat{s,ii}(:,1)>thr1,t_score_mat{s,ii}(:,2)<thr);
                % negpop = and(t_score_mat{s,ii}(:,2)>thr2,t_score_mat{s,ii}(:,1)<thr);
                
                thr = 0.075;
                thr1 = 0.075;
                thr2 = 0.075;
                pospop = and(w_score_mat{s,ii}(:,1)>thr1,w_score_mat{s,ii}(:,2)<thr);
                negpop = and(w_score_mat{s,ii}(:,2)>thr2,w_score_mat{s,ii}(:,1)<thr);

                % pospop = t_score_mat{s,ii}(:,1)>thr;
                % negpop = t_score_mat{s,ii}(:,2)>thr;

                vpop(s,ii,1) = sum(pospop);
                vpop(s,ii,2) = sum(negpop);

                neural_val_pos = neural_val(labels_val>binzc,pospop); % Pos pop coding pos
                neural_val_neg = neural_val(labels_val<binzc,negpop); % Neg pop coding neg
                neural_val_pos_cr = neural_val(labels_val>binzc,negpop); % Neg pop coding pos
                neural_val_neg_cr = neural_val(labels_val<binzc,pospop); % Pos pop coding neg

                labels_val_pos = labels_val(labels_val>binzc);
                labels_val_neg = labels_val(labels_val<binzc);

                % mdl = svmtrain(labels_val_neg, neural_val_neg,   '-s 3 -t 2 -c 1 -q');
                %  % Normalize training data
                % %  % Main/Mutual
                % [ rsa_P1(s,ii,1),~,rsa_P1t(s,ii,1)] = ARC_regress_wrapper(neural_val_neg,labels_val_neg,neural_val_pos_cr,labels_val_pos, 2,1);
                % [ rsa_P1(s,ii,2),~,rsa_P1t(s,ii,2)] = ARC_regress_wrapper(neural_val_pos,labels_val_pos,neural_val_neg_cr,labels_val_neg, 2,1);

                % Negpop
                [ rsa_P1(s,ii,1),~,rsa_P1t(s,ii,1)] = ARC_regress_wrapper(neural_val_pos_cr,labels_val_pos,neural_val_neg,labels_val_neg, 2,1);
                [ rsa_P1(s,ii,2),~,rsa_P1t(s,ii,2)] = ARC_regress_wrapper(neural_val_neg,labels_val_neg,neural_val_pos_cr,labels_val_pos, 2,1);
                % %
                % % Pospop
                error('check pop dist') % This maybe incorrect
                % [ rsa_P1(s,ii,1),~,rsa_P1t(s,ii,1)] = ARC_regress_wrapper(neural_val_neg,labels_val_neg,neural_val_pos_cr,labels_val_pos, 2,1);
                % [ rsa_P1(s,ii,2),~,rsa_P1t(s,ii,2)] = ARC_regress_wrapper(neural_val_pos_cr,labels_val_pos,neural_val_neg,labels_val_neg, 2,1);

                dfs(s,ii,1) = length(labels_val_pos);
                dfs(s,ii,2) = length(labels_val_neg);
                legender = {'Val--> Val+ in N-','Val+ -> Val- in N+'};

            otherwise
                error('Decoding type unspecified...')
        end
    end
end

%% Figure

% for plotsw = [true false]
%     if plotsw; rsa_plt = rsa_P1t; else; rsa_plt = rsa_P1; end
%     mkdir(savepath)
%     S_mat = squeeze(mean(rsa_plt));
%     S_err = squeeze(std(rsa_plt))./sqrt(3);
%     figure('Position',[0.5 0.5 640 480])
%     hold on
%     ngroups = size(S_mat, 1);
%     nbars = size(S_mat, 2);
%     b = bar(S_mat);
%     % b(1).FaceColor = [0 0.2470 0.9410];
%     % b(2).FaceColor = [0.3010 0.7450 0.9330];
%     % Calculating the width for each bar group
%     groupwidth = min(0.8, nbars/(nbars + 1.5));
%     x_m = [];
%     for i = 1:nbars
%         x = (1:ngroups) - groupwidth/2 + (2*i-1)*groupwidth/(2*nbars);
%         errorbar(x, S_mat(:,i), S_err(:,i), 'k.');
%         x_m = [x_m; x];
%     end
%     xticks(1:nanat)
%     xticklabels(anat_names);
%     c_s = {'r','g','b'}; % Data dots for subjects
%     for ii = 1:nanat % For bars for perceptual, chemical and combinations
%         for jj = 1:3
%             plot(x_m(:,ii),squeeze(rsa_plt(jj,ii,:)),c_s{jj})
%         end
%     end
%     if  plotsw
%         ylabel('Decoding accuracy (t)')
%     else
%         ylabel('Decoding accuracy (r)')
%     end
%     legend(legender)
%     if  plotsw
%         yline(1.65,'HandleVisibility','off')
%         yline(-1.65,'HandleVisibility','off')
%     end
%     % yline(r2t(0.025,4320),'HandleVisibility','off')
%     % yline(-r2t(0.025,4320),'HandleVisibility','off')
%     % legend({'Valence','Salience'})
%     clear modelmd_ modelmd modelmd2 S1_omat_vals S2_omat_vals unity M_anat M_sal M_val
%     p_values = arrayfun(@(r_eff) ARC_computePValueOneTailed(r_eff, binz, length(behav_ratings_)), rsa_P1);
%     if plotsw
%         savefig(fullfile(savepath,'ARC_decodingt'))
%         print(gcf,'-vector','-dsvg',[fullfile(savepath,'ARC_decodingt'),'.svg']) % svg
%     else
%         savefig(fullfile(savepath,'ARC_decoding'))
%         print(gcf,'-vector','-dsvg',[fullfile(savepath,'ARC_decoding'),'.svg']) % svg
%     end
% end

mkdir(savepath)
% savepath = pwd;
rsa_copy = rsa_P1;
p_values_3dt = ARC_decoding_pvals(rsa_P1, dfs);
% rsa_P1(vpop<10) = nan;
ARC_barplot_sig(rsa_P1,p_values_3dt)
xticks(1:nanat)
xticklabels(anat_names);
ylim([-0.2 0.1])
legend({'Val+ -> Val-','Val- -> Val+'})
% legend(legender)
savefig(fullfile(savepath,'ARC_decoding'))
print(gcf,'-vector','-dsvg',[fullfile(savepath,'ARC_decoding'),'.svg']) % svg
% print(gcf,'-vector','-dsvg',[fullfile(savepath,'ARC_decoding'),'.svg']) % svg
% p_values_3d2 = ARC_combinePValues3D(rsa_P1t, rsa_P1, dfs)

% p_values_3dz = ARC_decoding_pvals_zsc(rsa_P1, dfs)
% savepath = pwd;

save(fullfile(savepath,'ARC_decoding'))
toc

%% backup
% % mdl = svmtrain(labels_val_pos, neural_val_pos,   '-s 3 -t 2 -c 1 -q');
% % Normalize training data
% mu = mean(neural_val_pos);
% sigma = std(neural_val_pos);
% neural_val_pos_normalized = (neural_val_pos - mu) ./ sigma;
% % Normalize test data using the same mean and std from training data
% neural_val_neg_cr_normalized = (neural_val_neg_cr - mu) ./ sigma;
% mdl = ARC_regress_nested_mdl(neural_val_pos_normalized,labels_val_pos,2);
% predictions = svmpredict(  labels_val_neg , neural_val_neg_cr_normalized,  mdl, ' -q ');
% rsa_P1(s,ii,2) = fastcorr(predictions,labels_val_neg);
% [~,rsa_P1t(s,ii,2)] = ARC_r2t(rsa_P1(s,ii,2),length(labels_val_neg));

