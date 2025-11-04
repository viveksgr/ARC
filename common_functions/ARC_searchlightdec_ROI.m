function [dec_val1,dec_val2,df] = ARC_searchlightdec_ROI(cfg,anatmask,neural,labels,subj,r)

% Searchlight configuration
rget = cfg.rget;
sz = size(logical(anatmask));
ref_vox = round(sz/2);
[MX, MY, MZ] = ndgrid(1:sz(1),1:sz(2),1:sz(3));
radii = sqrt((MX-ref_vox(1)).^2 + (MY-ref_vox(2)).^2 + (MZ-ref_vox(3)).^2);
% prototype sphere index for radii<rget that can be used everywhere in the brain
radius_index = find(radii<rget) - sub2ind(sz,ref_vox(1),ref_vox(2),ref_vox(3)); % Neighbors index

% Gray mask configuration
lin_index = find(logical(anatmask));
nvox=length(lin_index);
% resultsvolume index
linindexconv = zeros(sz);
linindexconv(lin_index) = 1:length(lin_index);
if cfg.disc; bcn = cfg.binCentre; else; bcn = 0; end

if cfg.includeneut
    labels_disc = selectLess4AndHalf4(labels,bcn);
    labels_pos = labels(~labels_disc);
    labels_neg = labels(labels_disc);
else
    labels_pos = labels(labels>bcn);
    labels_neg = labels(labels<bcn);
end

dec_val1 = zeros(nvox,1);
dec_val2 = zeros(nvox,1);

% --- NEW: collect actual searchlight sizes within ROI ---
SLsize = zeros(nvox,1);
for cnt2 = 1:nvox
    indexindex2 = radius_index + lin_index(cnt2);
    indexindex2 = intersect(lin_index, indexindex2);
    indexindex2 = linindexconv(indexindex2);

    SLsize(cnt2) = numel(indexindex2);  % <<--- NEW: effective SL size at this center
    neural_val = double(neural(:,indexindex2));

    if (size(neural_val,1)>1)
        if cfg.includeneut
            neural_val_pos = neural_val(~labels_disc ,:);
            neural_val_neg = neural_val( labels_disc ,:);
        else
            neural_val_pos = neural_val(labels>bcn,:);
            neural_val_neg = neural_val(labels<bcn,:);
        end

        [dec_val1(cnt2,1)] = ARC_regress_wrapper(neural_val_pos,labels_pos, neural_val_neg, labels_neg, 4,1,false);
        [dec_val2(cnt2,1)] = ARC_regress_wrapper(neural_val_neg, labels_neg,neural_val_pos,labels_pos, 4,1,false);
    else
        dec_val1(cnt2,1) = 0;
        dec_val2(cnt2,1) = 0;
    end
end

df(1) = length(labels_pos);
df(2) = length(labels_neg);

validSL = SLsize > 0;
if any(validSL)
    meanSL = mean(SLsize(validSL));           % average #voxels per sphere inside ROI
    N_eff  = max(1, round(nvox / max(meanSL,1)));
else
    N_eff = max(1, nvox);                     % fallback (tiny ROI / degenerate)
end


subjDir  = fullfile(cfg.root,sprintf('ARC%02d',subj),'single');
roiName = cfg.ROIs.names{r};

if cfg.mapper
    dec_val1_map = unmasker(dec_val1,logical(anatmask));
    dec_val2_map = unmasker(dec_val2,logical(anatmask));
    write_reshaped_nifty(dec_val1_map, cfg.saveDir, false, fullfile(subjDir,cfg.maskFile), sprintf('ARC%02d_val1_RSAx1%s_pos',subj,roiName));
    write_reshaped_nifty(dec_val2_map, cfg.saveDir, false, fullfile(subjDir,cfg.maskFile), sprintf('ARC%02d_val2_RSAx1%s_pos',subj,roiName));
    write_reshaped_nifty(-dec_val1_map, cfg.saveDir, false, fullfile(subjDir,cfg.maskFile), sprintf('ARC%02d_val1_RSAx1%s_neg',subj,roiName));
    write_reshaped_nifty(-dec_val2_map, cfg.saveDir, false, fullfile(subjDir,cfg.maskFile), sprintf('ARC%02d_val2_RSAx1%s_neg',subj,roiName));
end