function matprod = ARC_matprod(ZI,ZJ)

ZI = ZI-nanmean(ZI(:)); 
ZJ = ZJ-nanmean(ZJ(:)); 

% % L2 normalize each column
% if min(min(ZI(:)),min(ZJ(:)))<0
%     error('Feed non negative vectors');
% end

ZI = ZI./sqrt(nansum(ZI.^2,'all'));
ZJ = ZJ./sqrt(nansum(ZJ.^2,'all'));
matprod  = (ZI.*ZJ);
