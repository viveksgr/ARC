function G = ARC_make_group_stats(wSub, tSub,cfg)
% Combine per-subject RSA outputs into group-level summaries

% Average p-values
nanat = length(cfg.anatNames);
rsa_Pavg = zeros(nanat,2);
rsa_P1wt = zeros(3,nanat,2);
for ii=1:nanat
    vec = cat(2,wSub{:,ii});
    vec_dist = cat(3,tSub{:,ii});
    for cond = 1:2    
        meaneff = nanmean(vec(cond+1,:));
        nulldist = squeeze(nanmean(vec_dist(:,cond+1,:),3));
        rsa_Pavg(ii,cond)= 1-invprctile(nulldist,meaneff)/100;
        rsa_P1wt(:,ii,cond) = vec(cond+1,:);
    end
end
G.rsa_Pavg = rsa_Pavg;

% Cross domain correlation
ARC_barplot_sig(rsa_P1wt,rsa_Pavg,true,~cfg.valenceSplit)
ylabel('RSA beta')
xticks(1:nanat)
xticklabels(cfg.anatNames)
if cfg.valenceSplit
    legend({'Val+','Val-'})
else
    legend({'Valence','Salience'})
end
% ylim(cfg.plotlims)
if cfg.verbose % save results figs etc
mkdir(cfg.saveRoot)
savefig(fullfile(cfg.saveRoot,'ARC_RSAwt'))
print(gcf,'-vector','-dsvg',[fullfile(cfg.saveRoot,'ARC_RSAwt'),'.svg']) % svg
print(fullfile(cfg.saveRoot,'ARC_RSAwt'),'-dpng')
end

end
