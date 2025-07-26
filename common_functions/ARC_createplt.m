function ARC_createplt(cfg,results_wt,results_dfs)
mkdir(cfg.saveDir);
beta = mean(results_wt,4,'omitmissing');
dfmat = mean(results_dfs,4,'omitmissing');
[pvals] = ARC_decoding_pvals(beta,dfmat);     % your existing helper
ARC_barplot_sig(beta,pvals)
xticks(1:cfg.nbROIs)
xticklabels(cfg.ROIs.names);
title(cfg.popChoice)
savefig(fullfile(cfg.saveDir,"ARC_decoding.fig"))
print(fullfile(cfg.saveDir,'ARC_decoding'),'-dpng')
print(gcf,'-vector','-dsvg',[fullfile(cfg.saveDir,'ARC_decoding'),'.svg']) % svg
save(fullfile(cfg.saveDir,'ARC_decoding.mat'));