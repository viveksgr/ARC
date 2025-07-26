function prctile_mat = ARC_prctile_mat_comp_prc(mat,cfg)

prctile_mat = cell(3,length(cfg.anatMasks));
for ss = 1:3
for ii = 1:length(cfg.anatMasks)

    % ARC_scatterDensity(w_mat{ii}(:,1),w_mat{ii}(:,2))
    cat1 = mat.pSub{ss, ii}.w_scores ;
    % cat2 = mat.pSub{ss, ii}.res_shuff.wt_mat;
    prctile_mat{ss,ii} = ARC_percentileAbsWeights( cat1 );

end
end