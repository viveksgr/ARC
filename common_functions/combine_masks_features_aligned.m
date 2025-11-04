function out = combine_masks_features_aligned(mask1, mask2, X1, X2, varargin)
% Combine two masks and features so that row order of out.X matches find(out.mask).
% Inputs:
%   mask1, mask2 : 3D logical/numeric arrays (same size) OR .nii filenames
%   X1 : [N1 x P] features for voxels where mask1==1  (row order must match find(mask1) or Index1)
%   X2 : [N2 x P] features for voxels where mask2==1  (row order must match find(mask2) or Index2)
%
% Name/Value:
%   'Merge'    : 'mean' (default) | 'first' | 'second' — policy for overlaps
%   'Index1'   : explicit linear indices for rows of X1 (default: find(mask1))
%   'Index2'   : explicit linear indices for rows of X2 (default: find(mask2))
%   'WriteUnion' : true/false (default false) — write union mask if inputs are NIfTI
%   'UnionFile'  : output filename (default next to mask1 as union_mask.nii)
%
% Output:
%   out.X        : [Nu x P] combined features (Nu = #voxels in union)
%   out.mask     : 3D logical union mask
%   out.idx      : [Nu x 1] linear indices (the order you get from find(out.mask))
%   out.has1/has2: logical vectors over out.idx (belongs to mask1/mask2)
%   out.map1/map2: row indices into X1/X2 for each voxel (0 if absent)
%   out.files    : paths used/written (if NIfTI)

ip = inputParser;
addParameter(ip,'Merge','mean',@(s) ismember(lower(s),{'mean','first','second'}));
addParameter(ip,'Index1',[]);
addParameter(ip,'Index2',[]);
addParameter(ip,'WriteUnion',false,@islogical);
addParameter(ip,'UnionFile','',@ischar);
addParameter(ip,'GMMask',[]);
parse(ip,varargin{:});
opt = ip.Results; merge = lower(opt.Merge);

% --- read masks (arrays or NIfTI) ---
isF1 = ischar(mask1) || isstring(mask1);
isF2 = ischar(mask2) || isstring(mask2);
if isF1, V1 = spm_vol(char(mask1)); m1 = spm_read_vols(V1)>0.5; else, m1 = mask1~=0; V1=[]; end
if isF2, V2 = spm_vol(char(mask2)); m2 = spm_read_vols(V2)>0.5; else, m2 = mask2~=0; V2=[]; end
if any(opt.GMMask(:)); m1 = and(m1,opt.GMMask); m2 = and(m2,opt.GMMask); end
if ~isequal(size(m1),size(m2)), error('Masks must have identical matrix size.'); end

% --- index lists used to map rows of X1/X2 ---
if isempty(opt.Index1), idx1 = find(m1); else, idx1 = opt.Index1(:); end
if isempty(opt.Index2), idx2 = find(m2); else, idx2 = opt.Index2(:); end
N1 = numel(idx1); N2 = numel(idx2);
if size(X1,1)~=N1, error('X1 rows (%d) must equal #indices in mask1 (%d).',size(X1,1),N1); end
if size(X2,1)~=N2, error('X2 rows (%d) must equal #indices in mask2 (%d).',size(X2,1),N2); end
P = size(X1,2);
if size(X2,2)~=P, error('X1 and X2 must have same #features.'); end

% --- union mask and its canonical order ---
uMask = m1 | m2;
uIdx  = find(uMask);                 % >>> THIS defines output row order <<<
Nu    = numel(uIdx);

% --- map union indices to X1/X2 row positions ---
[has1, loc1] = ismember(uIdx, idx1); % loc1: row in X1 (0 if absent)
[has2, loc2] = ismember(uIdx, idx2);

% --- build combined X aligned to uIdx order ---
X = nan(Nu, P);
if any(has1 & ~has2), X(has1 & ~has2, :) = X1(loc1(has1 & ~has2), :); end
if any(~has1 & has2), X(~has1 & has2, :) = X2(loc2(~has1 & has2), :); end
both = has1 & has2;
if any(both)
    switch merge
        case 'mean'
            X(both,:) = 0.5*(X1(loc1(both),:) + X2(loc2(both),:));
        case 'first'
            X(both,:) = X1(loc1(both),:);
        case 'second'
            X(both,:) = X2(loc2(both),:);
    end
end

% --- optionally write union mask ---
files = struct('mask1','', 'mask2','', 'union','');
if isF1, files.mask1 = char(mask1); end
if isF2, files.mask2 = char(mask2); end
if opt.WriteUnion && ~isempty(V1)
    Vout = V1; Vout.fname = opt.UnionFile;
    if isempty(Vout.fname)
        [p,~,~] = fileparts(files.mask1);
        Vout.fname = fullfile(p,'union_mask.nii');
    end
    Vout.descrip = 'union(mask1|mask2)';
    spm_write_vol(Vout, double(uMask));
    files.union = Vout.fname;
end

% --- pack ---
out = struct();
out.X     = X;            % rows aligned to out.idx = find(out.mask)
out.mask  = uMask;
out.idx   = uIdx(:);
out.has1  = has1;         out.has2 = has2;
out.map1  = loc1;         out.map2 = loc2;
out.files = files;
end

% Script stuff
% fpath = 'C:\Work\ARC\ARC\ARC02\single';
% mask1 = fullfile(fpath,'rwofc.nii');
% mask2 = fullfile(fpath,'rwvmpfc.nii');
% m3 = spm_read_vols(spm_vol(fullfile(fpath,'ARC3_anatgw.nii')));
% X1 = load(fullfile(fpath,'OFC','TYPED_FITHRF_GLMDENOISE_RR.mat'),'modelmd');
% X1 = squeeze(X1.modelmd);
% X2 = load(fullfile(fpath,'VMPFC','TYPED_FITHRF_GLMDENOISE_RR.mat'),'modelmd');
% X2 = squeeze(X2.modelmd);
% out = combine_masks_features_aligned(mask1, mask2, X1, X2, 'Merge','mean','GMMask',m3);
% savepath = fullfile(fpath,'ar_exmats');
% mkdir(savepath)
% save(fullfile(savepath,'ar_ex_mats.mat'),'out')
