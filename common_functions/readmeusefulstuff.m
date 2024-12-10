function readmeusefulstuff(fname)

% load(fullfile(pwd,fname),'rsa_Pavg','rsa_P1wt')
% round(rsa_Pavg,3)
% round(squeeze(mean(rsa_P1wt,1)),3)
load(fullfile(pwd,fname),'p_values_3dt','rsa_P1')
round(p_values_3dt,3)
round(squeeze(mean(rsa_P1,1)),3)
