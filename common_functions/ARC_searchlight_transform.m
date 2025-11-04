function [beta_mat,pvals] = ARC_searchlight_transform(beta,dfmat,pmat,mode)

nS = size(beta,1); % num subjects
nAnat = size(beta,2); % num ROIs

beta_mat = zeros(nS,nAnat,2);
for ss = 1:nS
    for ii = 1:nAnat
       
        % if and(ss==3,ii==4)
        %     'beep'
        % end
        if nargin<3
            [p1,~,n1,~] = ARC_bidirectionalFDR( beta{ss,ii,1}, dfmat(ss,ii,1));
            [p2,~,n2,~] = ARC_bidirectionalFDR( beta{ss,ii,2}, dfmat(ss,ii,2));
        else
            if strcmp(mode,'rmode')
                % fprintf('rmode')
                [p1,n1] = ARC_mean_posneg_thresholds( beta{ss,ii,1}, pmat(ss,1:2));
                [p2,n2] = ARC_mean_posneg_thresholds( beta{ss,ii,2}, pmat(ss,3:4));
                % [p2,~,n2,~] = ARC_mean_posneg_thresholds( beta{ss,ii,2}, dfmat(ss,ii,2),pmat(ss,3:4));

            elseif strcmp(mode,'pmode')
                 % fprintf('pmode')
                [p1,~,n1,~] = ARC_bidirectionalFDR( beta{ss,ii,1}, dfmat(ss,ii,1),pmat(ss,1:2));
                [p2,~,n2,~] = ARC_bidirectionalFDR( beta{ss,ii,2}, dfmat(ss,ii,2),pmat(ss,3:4));
            end
        end
     
        beta_mat(ss,ii,1) = (p1+p2)/2;
        beta_mat(ss,ii,2) = (n1+n2)/2;
             % beta_mat(ss,ii,1) = (p2);
        % beta_mat(ss,ii,2) = (n2);
    end
end

dls = ceil(cellfun(@(x) length(x),beta));

[p_vals, info] = binop_pvals_by_tail(beta_mat, squeeze(dls(:,:,1)), ...
    'AlphaVoxel', 0.025, 'TwoSidedVoxel', false, 'Alternative', 'greater');




pvals = squeeze(sum(beta_mat>4));
pvals(pvals==1)=0.05;
pvals(pvals==2)=0.01;
pvals(pvals==3)=0.001;
pvals(pvals==0)=0.5;

beta_mat = abs(beta_mat);