function ARC_create_behav_corr(behav)
behav_labs = {'Intensity';' Fishy';' Burnt';' Sour';' Decayed';...
    ' Musky';' Fruity';' Sweaty';' Cool';' Chemical';'Floral';' Sweet';' Warm';...
    ' Bakery';' Garlic';' Spicy';' Acidic';'Ammonia';'Edible';'Familiar'};

idx(1,:) = [1:18 19 19 19];
idx(2:3,:) = repmat([1:10 19 11:14 19 15 19 16:18],2,1);

ndisc = 20;
nodor = 160;
nS = 3;

behav_corrs = zeros(3,ndisc);
for ss = 1:nS 
    behav_rat_temp = [behav.behav(ss).ratings nan(nodor ,1)];
    behav_rat_sorted =  behav_rat_temp(:,idx(ss,:));
    kk = 0;
    for zz = [1 3:ndisc+1]
        kk = kk+1;
        behav_corrs(ss,kk) = fastcorr( behav_rat_sorted(:,2), behav_rat_sorted(:,zz) ); 
    end
end

figure()
behav_mean = nanmean(behav_corrs);
[behav_mean,argsort] = sort(behav_mean,'descend');
behav_corrs = behav_corrs(:,argsort);

behav_labs_sort = behav_labs(argsort);
bar(1:ndisc,behav_mean)
hold on
errorbar(1:ndisc,behav_mean,1.96*nanstd(behav_corrs)/sqrt(3),'.')
c = {'r.','g.','b.'};
for ss = 1:3
    plot(1:ndisc,behav_corrs(ss,:),c{ss})
end
xticks(1:ndisc)
xticklabels(behav_labs_sort)
xtickangle(90)
yline(r2t(0.025/18,160))
yline(-r2t(0.025/18,160))
savefig('pls_corr.fig')