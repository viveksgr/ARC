function out = invnorm_reslice_masks(maskDir, maskFiles, anatDir, varargin)
% invnorm_reslice_masks
%   Inverse-normalize atlas/mask(s) to subject native space, then reslice to a
%   chosen reference image (e.g., anat or first functional).
%
% Inputs
%   maskDir   : folder containing mask(s) defined in normalized/MNI space
%   maskFiles : char or cellstr of mask filenames (e.g., 'insula.nii' or {'insula.nii','DLPFC.nii'})
%   anatDir   : subject's anatomical folder (contains y_*.nii, anatomical .nii)
%
% Name/Value options
%   'AnatFilePattern' : cellstr patterns to find native anat (default {'s*.nii','means*.nii'})
%   'DefFilePattern'  : cellstr patterns to find forward def field y_*.nii (default {'y_*.nii','y_s*.nii','y_means*.nii'})
%   'InvDefName'      : inverse def filename to write (default 'y_inverse.nii')
%   'DoMakeInverse'   : true/false, make inverse if missing (default true)
%   'RefImage'        : full path to image used for reslicing (default: auto-pick in anatDir)
%   'RefPattern'      : patterns to find ref image if RefImage not given (default {'rf*.nii','s*.nii'})
%   'BB'              : bounding box for write (default [-78 -112 -85; 78 95 85])
%   'Vox'             : voxel size for write (default [1 1 1])
%   'InterpWrite'     : interpolation for write (0=nearest) (default 0; keeps masks binary-ish)
%   'InterpReslice'   : interpolation for reslice (default 0)
%   'CopyWToAnat'     : move w* masks to anatDir after write (default false)
%   'ReslicePrefix'   : prefix for resliced output (default 'r')
%
% Output (struct)
%   out.invDef   : path to inverse deformation
%   out.wMasks   : cellstr paths to inverse-normalized masks (w*)
%   out.rwMasks  : cellstr paths to resliced masks (r w*)
%   out.refImage : path used for reslicing
%
% Requires SPM12 on MATLAB path.

% --- defaults & parsing ---
ip = inputParser;
addParameter(ip,'AnatFilePattern', {'s*.nii','means*.nii'});
addParameter(ip,'DefFilePattern',  {'y_*.nii','y_s*.nii','y_means*.nii'});
addParameter(ip,'InvDefName',      'y_inverse.nii');
addParameter(ip,'DoMakeInverse',   true);
addParameter(ip,'RefImage',        '');
addParameter(ip,'RefPattern',      {'rf*.nii','s*.nii'});
addParameter(ip,'BB',              [-78 -112 -85; 78 95 85]);
addParameter(ip,'Vox',             [1 1 1]);
addParameter(ip,'InterpWrite',     0);
addParameter(ip,'InterpReslice',   0);
addParameter(ip,'CopyWToAnat',     false);
addParameter(ip,'ReslicePrefix',   'r');
parse(ip, varargin{:});
opt = ip.Results;

if ischar(maskFiles), maskFiles = {maskFiles}; end

% --- init SPM ---
spm('defaults','fmri');
try, spm_jobman('initcfg'); end

% --- locate anat file & forward deformation ---
anatFile = find_first(anatDir, opt.AnatFilePattern, true);
defFile  = find_first(anatDir, opt.DefFilePattern,  true);

% --- inverse deformation: create if missing ---
invDef = fullfile(anatDir, opt.InvDefName);
if ~exist(invDef,'file')
    if ~opt.DoMakeInverse
        error('Inverse deformation not found: %s', invDef);
    end
    matlabbatch = [];
    matlabbatch{1}.spm.util.defs.comp{1}.inv.comp{1}.def = {defFile};
    matlabbatch{1}.spm.util.defs.comp{1}.inv.space = {anatFile};
    matlabbatch{1}.spm.util.defs.out{1}.savedef.ofname = opt.InvDefName;
    matlabbatch{1}.spm.util.defs.out{1}.savedef.savedir.saveusr = {anatDir};
    spm_jobman('run', matlabbatch);
end

% --- apply inverse deformation to masks (write) ---
resampleList = cellfun(@(f) fullfile(maskDir,f), maskFiles, 'uni', false);
matlabbatch = [];
matlabbatch{1}.spm.spatial.normalise.write.subj.def      = {invDef};
matlabbatch{1}.spm.spatial.normalise.write.subj.resample = resampleList(:);
matlabbatch{1}.spm.spatial.normalise.write.woptions.bb   = opt.BB;
matlabbatch{1}.spm.spatial.normalise.write.woptions.vox  = opt.Vox;
matlabbatch{1}.spm.spatial.normalise.write.woptions.interp = opt.InterpWrite;
spm_jobman('run', matlabbatch);

% paths to written masks (prefix 'w')
wMasks = cellfun(@(f) fullfile(maskDir, ['w' strip_gz(f)]), maskFiles, 'uni', false);

% optionally move w* to anatDir
if opt.CopyWToAnat
    for i=1:numel(wMasks)
        tgt = fullfile(anatDir, get_filename(wMasks{i}));
        if ~strcmp(wMasks{i}, tgt)
            if exist(tgt,'file'), delete(tgt); end
            movefile(wMasks{i}, tgt);
            wMasks{i} = tgt;
        end
    end
end

% --- pick reslice reference ---
if ~isempty(opt.RefImage)
    refImage = opt.RefImage;
else
    refImage = find_first(anatDir, opt.RefPattern, true);
end
refCell = {sprintf('%s,1', refImage)};

% --- reslice to reference (coreg write) ---
matlabbatch = [];
matlabbatch{1}.spm.spatial.coreg.write.ref                = refCell;
matlabbatch{1}.spm.spatial.coreg.write.source             = wMasks(:);
matlabbatch{1}.spm.spatial.coreg.write.roptions.interp    = opt.InterpReslice;
matlabbatch{1}.spm.spatial.coreg.write.roptions.wrap      = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.write.roptions.mask      = 0;
matlabbatch{1}.spm.spatial.coreg.write.roptions.prefix    = opt.ReslicePrefix;
spm_jobman('run', matlabbatch);

% resliced outputs have prefix 'r' (default), applied to the *current* mask paths
rwMasks = prepend_prefix(wMasks, opt.ReslicePrefix);

% --- (optional) quick orientation sanity check ---
% v = spm_vol([{refImage}; rwMasks(:)]);
% spm_check_orientations([v{:}]);

% --- out ---
out = struct('invDef',invDef, 'wMasks',{wMasks}, 'rwMasks',{rwMasks}, 'refImage',refImage);

end

% ===== helpers =====
function f = find_first(dirpath, patterns, mustExist)
    if ischar(patterns), patterns = {patterns}; end
    f = '';
    for i=1:numel(patterns)
        dd = dir(fullfile(dirpath, patterns{i}));
        if ~isempty(dd)
            f = fullfile(dirpath, dd(1).name);
            break;
        end
    end
    if mustExist && isempty(f)
        error('File not found in %s for patterns: %s', dirpath, strjoin(patterns, ', '));
    end
end

function s = strip_gz(fname)
    % remove trailing .gz if present
    [p,n,e] = fileparts(fname);
    if strcmpi(e,'.gz')
        [~,n2,e2] = fileparts(n);
        s = [n2 e2];
    else
        s = [n e];
    end
end

function n = get_filename(pth)
    [~,n,ext] = fileparts(pth);
    n = [n ext];
end

function out = prepend_prefix(paths, prefix)
    out = cell(size(paths));
    for i=1:numel(paths)
        [p,n,e] = fileparts(paths{i});
        out{i} = fullfile(p, [prefix n e]);
    end
end
