function create_ARCscatter(w_score_mat,cfg)
% Scatter density
nanat = length((cfg.anatNames));
figure('Position', [0.5 0.5 1420 240])
% clims = {[0 5000],[0 6000],[0 60000],[0 30000]};
hold on
kk = 0;
for ii = 1:nanat
    kk = kk+1;
    subplot(1,nanat,kk)
    hold on
    % ARC_scatterDensity(w_mat{ii}(:,1),w_mat{ii}(:,2))

    w_mat = vertcat(w_score_mat{:,ii});
    w_mat = reshape(w_mat,[],2);
    histogram2( (w_mat(:,1)), (w_mat(:,2)),'DisplayStyle','tile','BinWidth',[range(w_mat(:,1))/40 range(w_mat(:,2))/40],'EdgeColor','none')

    % histogram2( w_mat{ii}(:,1), w_mat{ii}(:,2),'DisplayStyle','tile','BinWidth',[0.05 0.05],'EdgeColor','none')
    if cfg.valenceSplit
    xlabel('Val+ RSA beta')
    ylabel('Val- RSA beta')
    else
            xlabel('Valence beta')
            ylabel('Salience RSA beta')
    end
    xline(0)
    yline(0)
    % clim([0 6])
    % xlim([-0.4 1])
    % xticks([-0.4 0.3 1])
    % ylim([-0.4 1])
    % yticks([-0.4 0.3 1])
    colorbar
    % clim(clims{ii})
end
if cfg.verbose
savefig(fullfile(cfg.saveRoot,'ARC_dens2'))
print(fullfile(cfg.saveRoot,'ARC_dens2'),'-dpng')
print(gcf,'-vector','-dsvg',[fullfile(cfg.saveRoot,'ARC_dens2'),'.svg']) % svg
end