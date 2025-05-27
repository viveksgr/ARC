% Create single trials of fMRI runs using GLM Single package

%% General settings
linux_config = false;
s = 1;
nvol = 856; % 856 for S1, 876 for S2, S3;
anat_mask = 'ARC3_anatgw.nii';
mname = 'wm';
anat_name = {'wm'};
anat_area = {'rwm_main.nii'};

sess_i = 1;
sess_f = 3;
TRdur  = 0;
set_i = 1;
set_f = 4;
nruns = {[1:4] [1:4] [1:4] [1:4]};
TR = 1.4;
stimdur = 2;
nruns_T = length(horzcat(nruns{:}))*(sess_f-sess_i+1);
nfiles = nvol*nruns_T;

% % % Options for GLMsingle0
% opt.wantlss = 0;
% opt.wantlibrary = 0;
% opt.wantglmdenoise = 0;
% opt.wantfracridge  = 0;

if linux_config
    addpath('/home/vsh3681/spm12')
    addpath('/home/vsh3681/Scripts/imaging_analysis/common_functions')
    addpath('/home/vsh3681/Scripts/libsvm-3.14/matlab')
    addpath('/home/vsh3681/Scripts/Toolboxes/GLMsingle-main/GLMsingle-main/matlab')
    addpath('/home/vsh3681/Scripts/Toolboxes/fracridge-master/fracridge-master/matlab')
    addpath('/home/vsh3681/Scripts/Toolboxes/GLMsingle-main/GLMsingle-main/matlab/utilities')
    statpath = fullfile('/projects/p31178/ARC',sprintf('ARC%02d',s),'single');
    maskpaths = statpath;
    if s==3
        s=4;
    end
    datapath = fullfile('/projects/p31178/Data/NEMO',sprintf('NEMO_%02d',s),'/imaging/nii');
else
    addpath('C:\Work\Tools\GLMsingle-main\matlab')
    addpath('C:\Work\Tools\fracridge-master\fracridge-master\matlab')
    addpath('C:\Work\Tools\\GLMsingle-main\matlab\utilities')
    statpath = fullfile('C:\Work\ARC\ARC',sprintf('ARC%02d',s),'single');
    maskpaths = statpath;
    datapath = fullfile('C:\Work\NEMO Extended\Imaging\',sprintf('NEMO_%02d',s),'\nii');
end

% Load files
filestruct = cell(nfiles,1);
ii = 0;
rr = 0;
noise_files = {};
for set_ = set_i: set_f
    for sess_ = sess_i:sess_f
        for run_ = nruns{set_}
            rr = rr+1;
            path_ = fullfile(datapath, sprintf('set_%02d', set_)...
                ,sprintf('sess_%02d', sess_), sprintf('run_%02d', run_));
            n = dir(fullfile(path_, sprintf('fNEMO*.nii')));
            if s==1
                ns = dir(fullfile(path_,sprintf('nusiance_regresssors_NEMO_%02d_set_%02d_sess_%02d_run_%02d.txt',s,set_,sess_,run_)));
            else
                ns = dir(fullfile(path_,sprintf('nusiance_regresssors_NEMO%02d_set_%02d_sess_%02d_run_%02d.txt',s,set_,sess_,run_)));
            end
            noise_files{rr} =  load(fullfile(path_,ns(1).name));
            for i=1:nvol
                ii=ii+1;
                filestruct{ii,1} = fullfile(path_, sprintf('sr%s',n(i).name));
            end
        end
    end
end

% Load nuisance regressors
filename = fullfile(statpath, 'nuisance*.txt');
n = dir(filename);
mp = cat(1,noise_files{:});
% mp = load(fullfile(statpath, n(1).name));
% mp(:,35:end)=[];
tic
for aa = 1:length(anat_name)
    % Gray matter masks
    anatmask = (spm_read_vols(spm_vol(fullfile(maskpaths, anat_mask))));
    anatmask(isnan(anatmask)) = 0;
    anatmask = logical(anatmask);
    marea = (spm_read_vols(spm_vol(fullfile(maskpaths, anat_area{aa}))));
    marea(isnan(marea))=0;
    marea(marea<0.01) = 0;
    marea(marea>1) = 1;
    marea = logical(marea);
    if ~strcmp(anat_name{aa},'wm')
        anatmask = and(anatmask,marea);
    else
        anatmask = marea;
    end
    % Get data from all voxels in the mask
    Res_V = spm_vol(filestruct);
    [~, XYZmm] = spm_read_vols(Res_V{1});
    XYZvx = round(Res_V{1}.mat\[XYZmm; ones(1,size(XYZmm,2))]);
    Res_Vol = cell2mat(Res_V);
    anatmask1D = anatmask(:);
    vxl = XYZvx(:,anatmask1D);
    voxel_act = spm_get_data(Res_Vol,vxl)';
    voxel_act = single(voxel_act);
    [r,c] = find(isnan(voxel_act));
    voxel_act(r,:) = [];
    t_axis = (0:1:nfiles)*TR;

    % Design matrix
    load(fullfile(statpath,sprintf('conditions_NEMO%02d.mat',s)),'onsets')
    onsets_vec = vertcat(onsets{:});
    nodors = length(onsets);

    for tt = 1:length(TRdur)
        subdir = fullfile(statpath,mname, anat_name{aa},sprintf('TR_%02d',tt));
        mkdir(subdir)
        design_m = zeros(nfiles,nodors);
        for ii = 1:nodors
            onsets_ind = onsets{ii}+TRdur(tt); % Onset times for this cue
            zz_list = zeros(length(onsets_ind),1);
            for zz=1:length(zz_list)
                zz_list(zz) = find_nearest(t_axis,onsets_ind(zz));
            end
            design_m(zz_list,ii)=1;
        end

        % Change to runwise blocks
        design_m2 = cell(1,nruns_T);
        voxel_act2 = cell(1,nruns_T);
        mp2 = cell(1,nruns_T);
        for zz = 1:nruns_T
            idx = (zz-1)*nvol+1:(zz)*nvol;
            design_m2{zz} = design_m(idx,:);
            voxel_act2{zz} = voxel_act(:,idx);
            temp = mp(idx,:);
            mp2{zz} = temp;%(:,any(temp)); % Include non-zero regressors
        end
        opt.extraregressors = mp2;
        try
            results = GLMestimatesingletrial(design_m2,voxel_act2,stimdur,TR,subdir,opt);
        catch
            save(fullfile(subdir,'error_rep.mat'));
        end

    end
end
toc
