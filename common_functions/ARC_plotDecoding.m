function ARC_plotDecoding(beta,tstat,dfmat,cfg)
% Simple bar-plot with FDR-corrected p-values

mkdir(cfg.saveDir);
[pvals] = ARC_decoding_pvals(beta,dfmat);     % your existing helper

figure('Position',[100 100 640 480]); hold on
m   = mean(beta); se = std(beta)./sqrt(cfg.nbSubjects);
b   = bar(m); groupwidth = min(0.8,2/3);
for i = 1:2
    x = (1:cfg.nbROIs) - groupwidth/2 + (2*i-1)*groupwidth/4;
    errorbar(x,m(:,i),se(:,i),'k.');
end
xticks(1:cfg.nbROIs); xticklabels(cfg.ROIs.names);
ylabel('Decoding accuracy (r)'); legend({'Val+','Val-'});
title(sprintf('Decoding: %s population',cfg.popChoice));
print(fullfile(cfg.saveDir,'ARC_decoding'),'-dpng');
end
