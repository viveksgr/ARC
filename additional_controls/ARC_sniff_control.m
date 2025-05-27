%% RSA - Trial-wise

root = 'C:\Work';

settings_.nodor = 160;
settings_.wind = 7500; % Number of samples
settings_.nsniffcomp = 31;
settings_.loadvec = [3 4 9:settings_.nsniffcomp];
settings_.featorflessnot = true;
settings_.chem = false;
if  settings_.chem; delfeat = 1; else; delfeat = 0; end
rsa_pvals = zeros(3,4);

dirs = {fullfile(root ,'\SFP\sfp_behav_s01_correct');
    fullfile(root ,'\SFP\sfp_behav_s02_correct');
    fullfile(root ,'\SFP\sfp_behav_s04_correct')};

dirs2 = {fullfile(root,'ARC\ARC\ARC01\single');
    fullfile(root,'ARC\ARC\ARC02\single');
    fullfile(root,'ARC\ARC\ARC03\single')};

behav = load(fullfile('C:\Work\ARC\ARC\ARC','NEMO_perceptual2.mat'));
savepath = pwd;

sniff_corr = {};
figure('Position',[0 0 1280 320])
hold on
for ss = 1:length(dirs)
    % Breathing
    fprintf('Subject: %02d\n',ss)
    if ss==3; s2 = 4; else; s2 = ss; end
    statpath = dirs{ss};
    anatdir = dirs2{ss};
    % savepath = dirs3{ss};
    mkdir(savepath)

    load(fullfile(statpath,'sfp_feats_main.mat'))
    load(fullfile(statpath,'task_struct_trialwise.mat'))
    Fless_mat = vertcat(fless_mat{:});
    Fless_mat_pruned = Fless_mat(:,1:settings_.wind);

    Feat_mat_pruned = vertcat(feat_mat{:});
    Feat_mat_pruned =  Feat_mat_pruned(:,[settings_.loadvec]) ;
    Feat_mat_pruned(isnan(Feat_mat_pruned))=0;
    Feat_mat_pruned = zscore(Feat_mat_pruned,1);
   
    % Behavior
    onsets = load(fullfile(anatdir,sprintf('conditions_NEMO%02d.mat',s2)),'onsets');
    onsets = onsets.onsets;
    group_vec = cell(settings_.nodor,1);
    unity = [];
    for ii2 = 1:settings_.nodor
        group_vec{ii2} = ii2*ones(length(onsets{ii2}),1);
        unity = blkdiag(unity,ones(length(onsets{ii2})));
    end
    group_vec = vertcat(group_vec{:});
    [~,argsort] = sort(vertcat(onsets{:}));
    group_vec = group_vec(argsort);
    unity = unity(argsort,argsort);
    utl_mask = logical(triu(ones(length(group_vec)),1)); % All possible odors
    behav_ratings = behav.behav(ss).ratings(group_vec,:);
        
    inh_vol = Feat_mat_pruned(:,5);
    % inh_dur = Feat_mat_pruned(:,7);
    pls = behav_ratings(:,2);
    
    % 
    % subplot(3,3,ss)
    % hold on
    % plot(pls,inh_vol,'.')
    % rval = fastcorr(pls,inh_vol);
    % pval = r2p(rval,length(pls));
    % xlabel('Pleasantness')
    % ylabel('Sniff Dur')
    % title(sprintf('Sub: %02d (trials), r: %.3f, p: %.3f',ss,rval,pval))

    subplot(1,3,ss)
    hold on
    pls2 = splitapply(@mean,pls,group_vec);
    inh_vol2 = splitapply(@mean,inh_vol,group_vec);
    plot(pls2,inh_vol2,'.')
    X = [ones(length(pls2),1) pls2];
    beta = regress(inh_vol2,X);
    yreg = X*beta;
    plot(pls2,yreg)
    rval = fastcorr(pls2,inh_vol2);
    pval = r2p(rval,length(pls2));
    xlabel('Pleasantness')
    ylabel('Sniff Vol')
    title(sprintf('Sub: %02d (odors), r: %.3f, p: %.3f',ss,rval,pval))

    % subplot(3,3,6+ss)
    % hold on
    % 
    % % bin groups:
    pls = normalize(pls,'medianiqr');
    ms1 = min(pls);
    ms2 = max(pls);
    md = 0;
    edges = linspace( ms1,  ms2, 8);
    [~,~,bin] = histcounts(pls, edges);
    % 
    % pls2 = splitapply(@mean,pls,bin);
    % inh_vol2 = splitapply(@mean,inh_vol,bin);
    % plot(linspace(-1,1,7),inh_vol2,'.')
    % rval = fastcorr(pls2,inh_vol2);
    % pval = 2*r2p(rval,length(pls2));
    % xlabel('Pleasantness')
    % ylabel('Sniff Dur')
    title(sprintf('Sub: %02d (bins), r: %.3f, p: %.3f',ss,rval,pval))

    sniff_bin =    splitapply(@mean,Feat_mat_pruned,bin); 
    sniff_corr{ss} = corrcoef(sniff_bin');
end

savefig(fullfile(savepath,'fless_map'))
print(fullfile(savepath,'fless_map'), '-dpng')
SFP_clearLargeVariables(10)
% clear Fmat_1_m behav_mat Fless_mat Fless_mat_pruned A2_corr behav_corr int_corr pls_corr sess_run set_run task_run unity M_anat M_reg n_reg M_unstack fless_mat fless_mat_unn modelmd_ modelmd S_omat_vals utl_mask utl_mask2
save(fullfile(savepath,'ARC_RSA_fless'))
