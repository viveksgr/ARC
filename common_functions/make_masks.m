% make inverse normalize subject - job
% fpath = 'C:\Work\Templates\NeubertCingulateOrbitoFrontalParcellation\NeubertCingulateOrbitoFrontalParcellation\CingulateOrbitoFrontal_thr50_summaryimage_2mm.nii';
fpath = 'C:\Work\Templates\kahnt';
% ms_names = dir(fullfile(fpath,'ar*.nii'));
ms_names = dir(fullfile(fpath,'k_clus*.nii'));

maskDir = fpath;
for zz  = 1:length(ms_names)
ms_name = ms_names(zz).name;
fname1 = sprintf('w%s',ms_name);
fname2 = sprintf('rw%s',ms_name);

% Subject 1
anatDir = 'F:\NEMO_oldrepo\Older_repo-II\NEMO\NEMO_01\imaging\nii\master_anat';
savepath = 'C:\Work\ARC\ARC\ARC01\single';
out1 = invnorm_reslice_masks(maskDir, ms_name , anatDir, ...
       'AnatFilePattern', {'means*.nii'}, ...
       'DefFilePattern',  {'y_means*.nii'}, ...
       'CopyWToAnat', true, ...
       'RefPattern', {'rf*.nii','means*.nii'});
movefile(fullfile(anatDir,fname1), fullfile(savepath,fname1))
movefile(fullfile(anatDir,fname2 ), fullfile(savepath,fname2 ))

% Subject 2
anatDir = 'F:\NEMO_oldrepo\Older_repo-II\NEMO\NEMO_02\imaging\nii\anat';
savepath = 'C:\Work\ARC\ARC\ARC02\single';
out2 = invnorm_reslice_masks(maskDir, ms_name, anatDir, ...
       'AnatFilePattern', {'sNEMO02.nii'}, ...
       'DefFilePattern',  {'y_sNEMO02.nii'}, ...
       'CopyWToAnat', true,...
       'RefPattern',      {'rf*.nii','sNEMO*.nii'});
movefile(fullfile(anatDir,fname1), fullfile(savepath,fname1))
movefile(fullfile(anatDir,fname2 ), fullfile(savepath,fname2 ))

% % Subject 4
anatDir = 'F:\NEMO_oldrepo\Older_repo-II\NEMO\NEMO_04\imaging\nii\anat';
savepath = 'C:\Work\ARC\ARC\ARC03\single';
out4 = invnorm_reslice_masks(maskDir, ms_name , anatDir, ...
       'AnatFilePattern', {'sNEMO04.nii'}, ...
       'DefFilePattern',  {'y_sNEMO04.nii'}, ...
        'CopyWToAnat', true,...
       'RefPattern',      {'rf*.nii','sNEMO*.nii'});
movefile(fullfile(anatDir,fname1), fullfile(savepath,fname1))
movefile(fullfile(anatDir,fname2 ), fullfile(savepath,fname2 ))
end