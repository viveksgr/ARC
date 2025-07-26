function D = ARC_load_subject_data(sid,cfg)
D.behav = load(cfg.behavFile,'behav');        % whole struct
D.behav = D.behav.behav(sid);                 % single subject slice
% D.subjectDir = fullfile(cfg.mainRoot,sprintf('ARC%02d',sid),'single');
D.subjectDir = fullfile(cfg.mainRoot,sprintf('ARC%02d',sid),'single');
onsets = load(fullfile(D.subjectDir,sprintf('conditions_NEMO%02d.mat',sid)),'onsets');
D.onsets = onsets.onsets;

% Find trial-wise structure
group_vec = cell(cfg.nodor,1);
for ii2 = 1:cfg.nodor
    group_vec{ii2} = ii2*ones(length(D.onsets{ii2}),1);
end
group_vec = vertcat(group_vec{:});
[~,argsort] = sort(vertcat(D.onsets{:}));
D.group_vec = group_vec(argsort);

end
