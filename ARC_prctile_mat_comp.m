function prctile_mat = ARC_prctile_mat_comp(mat,cfg)

prctile_mat = cell(3,length(cfg.anatMasks));
for ss = 1:3
for ii = 1:length(cfg.anatMasks)

    % ARC_scatterDensity(w_mat{ii}(:,1),w_mat{ii}(:,2))
    cat1 = mat.pSub{ss, ii}.wt_mat ;
    cat2 = mat.pSub{ss, ii}.res_shuff.wt_mat;
    prctile_mat{ss,ii} = ARC_inversePercentiles( cat1, cat2 );

end
end