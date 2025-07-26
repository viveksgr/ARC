function [Ycell, model, partition, condition] = ...
           ARC_pcm_cell(metaRow, varargin)
% ARC_buildPCMinputs  Assemble PCM inputs for a single ROI
%
%   [Ycell, model, partition, condition] = ARC_buildPCMinputs(metaRow)
%
% INPUT
%   metaRow   1 × Nsubj struct array (all subjects for one ROI).
%             Each element must contain fields:
%                 .M     :  P × N_trials   (voxel × trial) or N × P
%                 .v1    :  T × 1 behavioural feature 1 (same length for all subj)
%                 .v2    :  T × 1 behavioural feature 2
%                 .bin   :  N_trials × 1 condition label (1..T)
%                 .run   :  N_trials × 1 run / partition label (1..R)
%
% VARARGIN    name-value pairs forwarded to pcm (e.g. 'scalePrior',…)
%
% OUTPUT
%   Ycell       {Nsubj} cell, each N_trials × P   (trials × voxels)
%   model       PCM feature model with 2 components (v1,v2)
%   partition   {Nsubj} cell, each N_trials × 1   (run index)
%   condition   {Nsubj} cell, each Z matrix N_trials × T
%
% EXAMPLE
%   % meta is Nsubj × Nroi
%   [Ycell,model,part,cond] = ARC_buildPCMinputs(meta(:,roiID));
%   [~,theta_hat] = pcm_fitModelGroup(Ycell,model,part,cond);
%
% Author: ChatGPT-o3   July 2025
% -------------------------------------------------------------------------

% ----- basic checks ------------------------------------------------------
nSubj = numel(metaRow);
assert(nSubj >= 1,'metaRow must be 1×Nsubj struct array');

% ----- build Y, Z, run vectors ------------------------------------------
Ycell     = cell(1,nSubj);
partition = cell(1,nSubj);
condition = cell(1,nSubj);

% first subject determines #conditions T
T = numel(unique(metaRow{1}.bin(:)));

for s = 1:nSubj
    S = metaRow{s};

    % Y:   trials × voxels
    Y = double(S.M');
    if size(Y,1) ~= numel(S.bin)
        error('Subject %d: #trials mismatch between M and bin.',s);
    end

    Ycell{s}       = Y;
    
    temp = datasample(1:2,length(S.run(:)));
    % partition{s}   = S.run(:);                               % N × 1
    partition{s}   = temp;  


    condIdx        = S.bin(:);                               % N × 1
    condIdx = condIdx(randperm(length(condIdx)));
    condition{s}   = pcm_indicatorMatrix('identity',condIdx);% N × T

    % sanity checks
    if numel(unique(condIdx)) ~= T
        error('Subject %d: condition set differs from first subject.',s);
    end
end

% ----- build two-feature model ------------------------------------------
v1 = zscore(metaRow{1}.v1(:));     % T × 1   ensure same across sub
% v1 = zscore(rand(7,1));
v2 = zscore(metaRow{1}.v2(:));

Ac(:,:,1) = [v1 , zeros(T,1)];     %   T × 2  (feature 1 in col-1)
Ac(:,:,2) = [zeros(T,1) , v2];     %   T × 2  (feature 2 in col-2)

model.type        = 'feature';
model.name        = 'R1+R2 PCM';
model.numGparams  = 2;
model.Ac          = Ac;            % T × 2 × 2
model.theta0      = [0 ; 0];       % log-scale start (optional)

% % ----- optional overrides via varargin -----------------------------------
% if ~isempty(varargin)
%     model = pcm_setparams(model,varargin{:});   % helper from toolbox
% end

[~,wt] = pcm_fitModelGroup(Ycell, model, partition, condition , ...
                      'runEffect','none');   % safest if you have ≥1 run
% [~,wt] = pcm_fitModelIndivid(Ycell, model, partition, condition , ...
                      % 'runEffect','none');   % safest if you have ≥1 run

end
