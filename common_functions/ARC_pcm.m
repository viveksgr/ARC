function [res] = ARC_pcm(v1,v2,M,bin,run,cfg)

if cfg.runpcmcollapse
    res.v1 = v1;
    res.v2 = v2;
    res.M = M;
    res.bin = bin;
    res.run = run;
else
Ysingle = double(M');            % N × P
Ycell    = {Ysingle};                        % wrap into cell for pcm_*

% ------------------------------------------------------------------------
% 2) Partition vector : which run / fold each trial belongs to
% ------------------------------------------------------------------------
runIdx   = run;             % N × 1   (integers 1..R)
partition = {runIdx};                        % cell

% ------------------------------------------------------------------------
% 3) Condition vector : which stimulus / condition for each trial
% ------------------------------------------------------------------------
condIdx  = bin;       % N × 1   (integers 1..T)
condition = {condIdx};                       % cell

% If you already have a design matrix Z (N×T) use condition = {Z};
% behavioural regressors (T × 1)
v1  = zscore(v1(:));          % pleasantness feature
v2  = zscore(v2(:));          % salience    feature
T   = numel(v1);              % #conditions
% ------------------------------------------------------------------
% stack the two single-column feature matrices into Ac(:,:,k)
% – trick:  put each feature in its *own column*;
%   initialise the other column with zeros, then let θ scale them
% if ~cfg.valenceSplit
Ac1            = [v1  zeros(T,1)];   % 1st θ scales col-1  (R1)
Ac2            = [zeros(T,1)  v2];   % 2nd θ scales col-2  (R2)
% Ac1            = [v1];   % 1st θ scales col-1  (R1)
% Ac2            = [v2];   % 2nd θ scales col-2  (R2)
model                 = struct;
model.type            = 'feature';
model.numGparams      = 2;           % θ1, θ2
model.Ac              = cat(3,Ac1,Ac2);   % T × 2 × 2
model.name            = 'Feature-PCM (R1 & R2)';

% M.type = ‘component’;
% M.numGparams = 2;
% M.Gc(:,:,1) = Model(1).G_cent;
% M.Gc(:,:,2) = Model(2).G_cent;
% M.name = 'muscle + usage';

% else
% % Z-score
% v1 = zscore(v1);
% v2 = zscore(v2);
% 
% % Remove shared variance symmetrically
% v2u = v2 - (v1'*v2)/(v1'*v1)*v1;   % residualise v2 wrt v1
% v1u = v1 - (v2'*v1)/(v2'*v2)*v2;   % optional symmetric step
% 
% % Build 3-feature model: v1u, v2u, vShare       
% vShare = (v1 + v2)/2;
% 
% Ac(:,:,1) = [v1u , 0    , 0];
% Ac(:,:,2) = [0   , v2u  , 0];
% Ac(:,:,3) = [0   , 0    , vShare];
% 
% model.type       = 'feature';
% model.numGparams = 3;
% model.Ac         = Ac;
% end

Ysingle = Ycell{1};
% theta0_model = log(0.5) * ones(model.numGparams,1);   % here numGparams = 2

% ------------------------------------------------------------------
% 2) ask the toolbox for a starting θ
% % ------------------------------------------------------------------
% [wtx,theta_hat,~] = ...
%     pcm_fitModelGroup(Ycell, model, partition, condition , ...
%                       'runEffect','fixed');   % safest if you have ≥1 run
% 
% 3. build a fresh design matrix Z_perm and hand it to PCM explicitly
Z_perm = pcm_indicatorMatrix('identity', bin);
% 4. supply the design matrix itself, not just the vector, e.g.
condition_perm = {Z_perm};    

% 
[wtx,theta_hat,~] = pcm_fitModelIndivid(Ycell, model, partition, condition_perm , ...
                      'runEffect','none');   % safest if you have ≥1 run

% partition{1} = datasample(1:2,length(partition{1}));
% [wtx,theta_hat,~] = pcm_fitModelIndivid(Ycell, model, partition, condition , ...
%                       'runEffect','none');   % safest if you have ≥1 run



Z = pcm_indicatorMatrix('identity', condIdx);   % N × T
X = [];        

[Wfeat,varW] = pcm_estimateW( model,           ...   % feature model !
                           theta_hat{1},    ...   % θ̂  (vector)
                           Ycell{1},        ...   % Y  (N×P)
                           Z, X,            ...   % design, fixed eff.
                           'runEffect',[]);       % matches runEffect
betaw = theta_hat{1};
voxelWeights = Wfeat.';   
res.w_scores = voxelWeights;
res.t_corr = fastcorr(voxelWeights(:,1),voxelWeights(:,2));
res.beta = betaw(1:model.numGparams);
end

% [A, varW] = pcm_getPosterior(model, theta_hat{1},  Ycell{1}, Z, X, 'runEffect',[]);
