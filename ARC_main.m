%--------------------------------------------------------------
%  main_RSA.m   — driver for Representational‐Similarity analysis
%--------------------------------------------------------------
clear; clc;

cfg = ARC_make_default_config();                % one struct with all toggles
results = ARC_run_RSA_full(cfg);                % do the work
save(fullfile(cfg.saveRoot,'RSA_Group.mat'),'-struct','results');
disp('✓  RSA finished and saved');
