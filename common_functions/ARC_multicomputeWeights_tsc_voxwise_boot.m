function [t_corr, t_scores, proportions, w_scores, t_corr_sig, wt_mat] = ARC_multicomputeWeights_tsc_voxwise_boot(m1, m2, anatmat,seq)
% m1, m2: bin x bin matrices
% anatmat: N_voxels x bin x B bootstraps

if nargin<4
    seq = false;
end

[nvox, binz, nboot] = size(anatmat);

if mod(binz,2)==0
    binzpart1 = binz/2;
    binzpart2 = binzpart1+1;
else
    binzpart1 = (binz+1)/2;
    binzpart2 = binzpart1;
end

t_scores = zeros(nvox, 2);  % Averaged across bootstraps
w_scores = zeros(nvox, 2);  % Averaged across bootstraps

% Utility mask (same as before)
utl_mask1 = logical(blkdiag(zeros(binz-binzpart1), ones(binzpart1)));
utl_mask_blk = flipud(utl_mask1);
utl_mask_blk(binzpart1, binzpart1) = false;

wt_mat = zeros(nvox,nboot, 2);
for ii = 1:nvox
    wt_all = zeros(nboot, 2);
    t_all = zeros(nboot, 2);
    for b = 1:nboot
        % modelmd_binned = zscore(squeeze(anatmat(ii,:,b)));
        modelmd_binned = squeeze(anatmat(ii,:,b));
        anatcorr = 1 - abs(modelmd_binned - modelmd_binned');

        % x = squeeze(anatmat(ii,:,b));
        % % Normalize to [0,1] before computing "1 − abs" similarity
        % x_min = min(x);
        % x_rng = max(x) - x_min;
        % if x_rng > eps
        %     x_norm = (x - x_min) / x_rng;
        % else
        %     x_norm = zeros(size(x));  % flat response
        % end
        % anatcorr = 1 - abs(x_norm - x_norm');
        
        % anatcorr = modelmd_binned'*modelmd_binned;
        % anatcorr = gaussdist(modelmd_binned);

        anatcorr = squeeze(anatcorr);
        if seq
            y = anatcorr(utl_mask_blk);
            x1 = m1(utl_mask_blk);
            x2 = m2(utl_mask_blk);
            y1 = regressmeout(y',x2')';
            y2 = regressmeout(y',x1')';
            [wt, t] = ARC_multicomputeWeights_tsc(x1, y1);
            wt_all(b, 1) = wt(2);
            t_all(b, 1) = t(2);
            [wt, t] = ARC_multicomputeWeights_tsc(x2, y2);
            wt_all(b, 2) = wt(2);
            t_all(b, 2) = t(2);       
        else
            [wt, t] = ARC_multicomputeWeights_tsc([m1(utl_mask_blk),m2(utl_mask_blk)], anatcorr(utl_mask_blk));
              % [wt, t] = ARC_multicomputeWeights_tsc([m1(utl_mask2),m2(utl_mask2)], anatcorr(utl_mask2));
           
            wt_all(b, :) = wt(2:end);
            t_all(b, :) = t(2:end);
        end
      
    end
    wt_mat(ii,:,:) = wt_all;
    w_scores(ii,:) = mean(wt_all, 1);

     t_scores(ii,:) = mean(t_all, 1);
end

t_corr = corr(w_scores(:,1), w_scores(:,2));
t_corr_sig = ARC_r2t(t_corr, nvox);
% proportions = computeProportions(t_scores, tinv(0.975, nvox));
proportions = computeProportions(t_scores, 0.3,0.075);
end

function proportions = computeProportions(t_scores, thr,thr2)
    % Calculate conditions
    col1Above = t_scores(:,1) > thr & t_scores(:,2) <= thr2; % Col 1 above thr, Col 2 not
    col2Above = t_scores(:,2) > thr & t_scores(:,1) <= thr2; % Col 2 above thr, Col 1 not
    bothAbove = t_scores(:,1) > thr & t_scores(:,2) > thr;  % Both cols above thr
    
    % Calculate proportions
    totalRows = size(t_scores, 1);
    propCol1Above = sum(col1Above) / totalRows;
    propCol2Above = sum(col2Above) / totalRows;
    propBothAbove = sum(bothAbove) / totalRows;
    
    % Combine proportions into a 1x3 vector
    proportions = [propCol1Above, propCol2Above, propBothAbove];
end

function S = gaussdist(vec)
nbins = length(vec);
vec = vec(:); % Ensure column vector
S = zeros(nbins);
sigma = 1;
for i = 1:nbins
    for j = 1:nbins
        if i ~= j
            dist = vec(i) - vec(j);
            S(i,j) = exp(-(dist^2) / (2 * sigma^2));
        end
    end
end
end

% % Debug plots
% edge = -0.3:0.025:0.7;
% figure()
% hold on
% subplot(1,2,1)
% hold on
% histogram(w_scores(:,1),edge)
% histogram(w_scores2(:,1),edge)
% subplot(1,2,2)
% hold on
% histogram(w_scores(:,2),edge)
% histogram(w_scores2(:,2),edge)