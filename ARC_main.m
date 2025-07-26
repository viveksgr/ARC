% %--------------------------------------------------------------
% %  main_RSA.m   — driver for Representational‐Similarity analysis
% %--------------------------------------------------------------
% % tic

 % clear; clc;
cfg = ARC_make_default_config();                % one struct with all toggles
% cfg.seed = 2;
results = ARC_run_RSA_full(cfg);                % do the work

% % % disp('✓  RSA finished and saved');
% % toc
% % 
% % tic
% % cfg = ARC_makeDecodingCfg('C:\Work\ARC\ARC');
% % ARC_runDecoding(cfg);
% % toc
% 
% %% Iterative RSAtic
% clear; clc;
% nseeds = 1;
% res_wts = cell(3,4,nseeds);
% res_ts = cell(3,4,nseeds);
% results_cell = cell(nseeds,1);
% tic
% for seed = 1:nseeds
%     fprintf('Seed: %02d\n',seed)
%     cfg = ARC_make_default_config();
%     cfg.seed = seed;
%     results_cell{seed} = ARC_run_RSA_full(cfg);
%     res_wts(:,:,seed) = results_cell{seed}.wSub;
%     res_ts(:,:,seed) = results_cell{seed}.tSub;
% end
% toc
% wSub = ARC_cellAvg3D(res_wts);
% tSub = ARC_cellAvg3D(res_ts);
% group = ARC_make_group_stats(wSub, tSub,cfg);
% mkdir(cfg.saveRoot)
% save(fullfile(cfg.saveRoot,'wt_mats.mat'))
% 
% 
% %% Iterative decoding
% tic
% nseeds = 20;
% results_wt = zeros(3,4,2,nseeds);
% results_dfs = zeros(3,4,2,nseeds);
% for seed = 1:nseeds
%   fprintf('Seed: %02d\n',seed)
% cfg = ARC_makeDecodingCfg('C:\Work\ARC\ARC',seed);
% tmp = ARC_runDecoding(cfg);
% results_wt(:,:,:,seed) = tmp.beta;
% results_dfs(:,:,:,seed) = tmp.dfmat;
% end
% toc
% 
% mkdir(cfg.saveDir);
% beta = mean(results_wt,4,'omitmissing');
% dfmat = mean(results_dfs,4,'omitmissing');
% [pvals] = ARC_decoding_pvals(beta,dfmat);     % your existing helper
% ARC_barplot_sig(beta,pvals)
% xticks(1:cfg.nbROIs)
% xticklabels(cfg.ROIs.names);
% title(cfg.popChoice)
% savefig(fullfile(cfg.saveDir,"ARC_decoding.fig"))
% print(fullfile(cfg.saveDir,'ARC_decoding'),'-dpng')
% print(gcf,'-vector','-dsvg',[fullfile(cfg.saveDir,'ARC_decoding'),'.svg']) % svg
% save(fullfile(cfg.saveDir,'ARC_decoding.mat'));
% 
%% RSA + decoding
% clear; clc;
nseeds = 20;

% results_cell = cell(nseeds,1);
tic
results_wt = cell(3,nseeds);
results_dfs = cell(3,nseeds);
rootDir = 'C:\Work\ARC\ARC';
for seed = 1:nseeds
    fprintf('Seed: %02d\n',seed)
    cfg = ARC_make_default_config();
    cfg.seed = seed;
    results_cell{seed} = ARC_run_RSA_full(cfg);
    % prctile_mat = ARC_prctile_mat_comp(results_cell{seed},cfg);
    prctile_mat = ARC_prctile_mat_comp_prc(results_cell{seed},cfg);
    
    % prctile_mat = 
    % cfg_pos = ARC_makeDecodingCfg(rootDir,seed);
    % cfg_pos.popChoice = 'pos';
    % cfg_pos.samepop = true;
    % cfg_pos.saveDir      = fullfile(rootDir,'Decoding','prcseed','samepop');
    % cfg_pos.prctile_mat = prctile_mat;
    % tmp = ARC_runDecoding(cfg_pos);
    % results_wt{1,seed} = tmp.beta;
    % results_dfs{1,seed} = tmp.dfmat;

    cfg_neg = ARC_makeDecodingCfg(rootDir,seed);
    cfg_neg.popChoice = 'neg';
    cfg_neg.samepop = false;                                
    cfg_neg.saveDir      = fullfile(rootDir,'Decoding','temp','crossneg');
    cfg_neg.prctile_mat = prctile_mat;
    tmp = ARC_runDecoding(cfg_neg);
    results_wt{2,seed} = tmp.beta;
    results_dfs{2,seed} = tmp.dfmat;

    cfg_mut = ARC_makeDecodingCfg(rootDir,seed);
    cfg_mut.popChoice = 'mut';                               
    cfg_mut.saveDir      = fullfile(rootDir,'Decoding','temp','crossmut_lda');
    cfg_mut.samepop = false;
    cfg_mut.prctile_mat = prctile_mat;
    tmp = ARC_runDecoding(cfg_mut);
    results_wt{3,seed} = tmp.beta;
    results_dfs{3,seed} = tmp.dfmat;
end
% 
results_wt_i = cat(4,results_wt{3,:});
results_dfs_i = cat(4,results_dfs{3,:});
ARC_createplt(cfg_mut,results_wt_i,results_dfs_i)
% % toc
% % wSub = ARC_cellAvg3D(res_wts);
% % tSub = ARC_cellAvg3D(res_ts);
% % group = ARC_make_group_stats(wSub, tSub,cfg);
% % mkdir(cfg.saveRoot)
% % save(fullfile(cfg.saveRoot,'wt_mats.mat'))
% 
% 
%% Cross-decoding searchlight
rootDir = 'C:\Work\ARC\ARC';
cfg_cross = ARC_makeDecodingCfg(rootDir);
cfg_cross.saveDir      = fullfile(rootDir,'Decoding','cross_searchl_2_temp');
cfg_cross.includeneut =false;
ARC_runDecoding_searchl(cfg_cross)