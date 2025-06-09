function [val1, val2, maskMat] = ARC_build_behav_RSM(cfg);

binz = cfg.nBin;
binzpart1 = cfg.binzPart1;
binzpart2 = cfg.binzPart2;

% Contruction of theoretical RSMs
val_sc = linspace(-1,1,binz);
sal_sc = abs(val_sc);
valp = val_sc;
valp(1:binzpart1)=0;
valn = val_sc;
valn(binzpart2:end)=0;

val_mat = 1-abs(val_sc-val_sc');
valp_mat = 1-abs(valp-valp');
valn_mat = 1-abs(valn-valn');
sal_mat = 1-abs(sal_sc-sal_sc');

if cfg.runSqDist
    valp_mat = 1-abs(valp-valp').^2;
    valn_mat = 1-abs(valn-valn').^2;
end

utl_mask1 = logical(blkdiag(zeros(binz-binzpart1),ones(binzpart1)));
utl_mask2 = logical(triu(ones(length(val_mat)),1)); % All possible odors
utl_mask = and(utl_mask1,utl_mask2);
utl_mask_blk = flipud(utl_mask1);
utl_mask_blk( binzpart1 , binzpart1 )=false;

if cfg.valenceSplit
    val1 = valp_mat;
    val2 = valn_mat;
    maskMat = utl_mask2;
else
    val1 = valp_mat;
    val2 = valn_mat;
    maskMat = utl_mask;
end
