function [beta_mat,pvals] = ARC_searchlight_transform(beta,dfmat);

nS = size(beta,1); % num subjects
nAnat = size(beta,2); % num ROIs

beta_mat = zeros(nS,nAnat,2);
for ss = 1:nS
    for ii = 1:nAnat
        cellanat = beta{ss,ii,1};
        df1 = dfmat(ss,ii,1);   
        [p1,~,n1,~] = ARC_bidirectionalFDR( beta{ss,ii,1}, dfmat(ss,ii,1));
        [p2,~,n2,~] = ARC_bidirectionalFDR( beta{ss,ii,2}, dfmat(ss,ii,2));
        beta_mat(ss,ii,1) = (p1+p2)/2;
        beta_mat(ss,ii,2) = (n1+n2)/2;
    end
end

pvals = squeeze(sum(beta_mat>4));
pvals(pvals==1)=0.05;
pvals(pvals==2)=0.01;
pvals(pvals==3)=0.001;
pvals(pvals==0)=0.5;