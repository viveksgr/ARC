%--------------------------------------------------------------
%  main_RSA.m   — driver for Representational‐Similarity analysis
%--------------------------------------------------------------
tic
clear; clc;
cfg = ARC_make_default_config();                % one struct with all toggles
results = ARC_run_RSA_full(cfg);                % do the work
save(fullfile(cfg.saveRoot,'ARC_RSA.mat'));
disp('✓  RSA finished and saved');
toc


tic
cfg = ARC_makeDecodingCfg('C:\Work\ARC\ARC');
ARC_runDecoding(cfg);
toc


cfg = ARC_makeDecodingCfg('C:\Work\ARC\ARC');
ARC_runDecoding(cfg);

