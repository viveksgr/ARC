% Read the file as a string
filename = 'C:\Work\ARC\ARC_DNN\sp_sal-12-3-9\variables.json'; % replace with the path to your file
fileText = fileread(filename);

salier = false;

% Decode the JSON string
dataStruct = jsondecode(fileText);

pmat = dataStruct.p_mat;
h2_layer = dataStruct.h2_mat;
nsamp = 5;

%% Positive correlations
figure('Position',[0 0 1220 670])
k = 0;
for ss = 1:3
    C = load(fullfile('C:\Work\ARC\ARC_DNN\inpfiles',sprintf("C%01d.txt",ss)));
    val = C(:,end);
    if salier; val = abs(val); end
    a_ind = squeeze(pmat(ss,3,:));
    [~,argsort] = sort(a_ind);
    [~,argsort2] = sort(abs(a_ind));
    for ii = 1:nsamp
        k =  k+1;
        subplot(3,5,k)
        hold on
        manif2 = squeeze(h2_layer(ss,argsort(end-ii+1),:,:));
        % manif2 = squeeze(h2_layer(ss,argsort2(ii),:,:));
        scatter3(manif2(:,1), manif2(:,2), manif2(:,3), 25, val,'filled'); % 'filled' will fill the markers with color
        xlabel('X-axis');
        ylabel('Y-axis');
        zlabel('Z-axis');
        grid on
        view([1 1 1])
    end
end
title('Max correlations')
savefig(sprintf('manif_max'))
print(sprintf('manif_max'),'-dpng')

%% Histograms of dims
figure('Position',[0 0 1220 670])
k = 0;
for ss = 1:3
    C = load(fullfile('C:\Work\ARC\ARC_DNN\inpfiles',sprintf("C%01d.txt",ss)));
    val = C(:,end);
    if salier; val = abs(val); end
    a_ind = squeeze(pmat(ss,3,:));
    [~,argsort] = sort(a_ind);
    [~,argsort2] = sort(abs(a_ind));
    for ii = 1:3
        k =  k+1;
        subplot(3,3,k)
        hold on
        manif2 = squeeze(h2_layer(ss,argsort(end-ii+1),:,:));
        % manif2 = squeeze(h2_layer(ss,argsort2(ii),:,:));
        scatter3(manif2(:,1), manif2(:,2), manif2(:,3), 25, val,'filled'); % 'filled' will fill the markers with color
        xlabel('X-axis');
        ylabel('Y-axis');
        zlabel('Z-axis');
        grid on
        view([1 1 1])
    end
end
title('Max correlations')
savefig(sprintf('manif_max'))
print(sprintf('manif_max'),'-dpng')





% %% NUll correlations
% figure('Position',[0 0 1220 670])
% k = 0;
% for ss = 1:3
%     C = load(fullfile('C:\Work\ARC\ARC_DNN\inpfiles',sprintf("C%01d.txt",ss)));
%     val = C(:,end);
%     if salier; val = abs(val); end
%     a_ind = squeeze(pmat(ss,3,:));
%     [~,argsort] = sort(a_ind);
%     [~,argsort2] = sort(abs(a_ind));
%     for ii = 1:nsamp
%         k =  k+1;
%         subplot(3,5,k)
%         hold on
%         % manif2 = squeeze(h2_layer(ss,argsort(end-ii+1),:,:));
%         manif2 = squeeze(h2_layer(ss,argsort2(ii),:,:));
%         scatter3(manif2(:,1), manif2(:,2), manif2(:,3), 25, val,'filled'); % 'filled' will fill the markers with color
%         xlabel('X-axis');
%         ylabel('Y-axis');
%         zlabel('Z-axis');
%         grid on
%         view([1 1 1])
%     end
% end
% title('Null correlations')
% savefig(sprintf('manif_null'))
% print(sprintf('manif_null'),'-dpng')
% 
% %% Negative correlations
% figure('Position',[0 0 1220 670])
% k = 0;
% for ss = 1:3
%     C = load(fullfile('C:\Work\ARC\ARC_DNN\inpfiles',sprintf("C%01d.txt",ss)));
%     val = C(:,end);
%     if salier; val = abs(val); end
%     a_ind = squeeze(pmat(ss,3,:));
%     [a_ind_sort,argsort] = sort(a_ind);
%     [~,argsort2] = sort(abs(a_ind));
%     for ii = 1:nsamp
%         k =  k+1;
%         subplot(3,5,k)
%         hold on
%         manif2 = squeeze(h2_layer(ss,argsort(ii),:,:));
%         % manif2 = squeeze(h2_layer(ss,argsort2(ii),:,:));
%         scatter3(manif2(:,1), manif2(:,2), manif2(:,3), 25, val,'filled'); % 'filled' will fill the markers with color
%         xlabel('X-axis');
%         ylabel('Y-axis');
%         zlabel('Z-axis');
%         grid on
%         view([1 1 1])
%     end
% end
% title('Min correlations')
% savefig(sprintf('manif_min'))
% print(sprintf('manif_min'),'-dpng')

%% Correlation of Maps

root = 'C:\Work\ARC\ARC';
dirs = {fullfile(root,'\ARC01\mediation');
    fullfile(root,'\ARC02\mediation');
    fullfile(root,'\ARC03\mediation')};

behavP = load(fullfile(root,'ARC','NEMO_perceptual2.mat'));
behavC = load(fullfile(root,'ARC','NEMO_chemical2.mat'));
percepter = true;
utl_mask = logical(triu(ones(160),1)); % All possible odors

DNNner = true;
DNNloc = 'C:\Work\ARC\ARC_DNN\sp_val-12-3-7';
load(fullfile(DNNloc,'h2.mat'))
DNNloc = 'C:\Work\ARC\ARC_DNN\sp_sal-12-3-9';
load(fullfile(DNNloc,'h2_sal.mat'))

anat_cell = {};
matr = zeros(3,2);
% Subject - index
for s = [1 2 3] % Subject
    fprintf('Subject: %02d\n',s)
    anatdir = fullfile(root,sprintf('ARC%02d',s),'single');
    % Construction of median split
    behav_ratingsC = behavC.behav(s).ratings(:,1:20); % take only first 20 components of chemical space
    behav_ratingsP = behavP.behav(s).ratings(:,2); % 
    if percepter; X_mat = behavP.behav(s).ratings(:,[1 3:end]); else; X_mat =  behav_ratingsC; end

    [~, ~, p_values_val] = bootstrapRidge(X_mat, behav_ratingsP, 1000, 0.1);
    [~, ~, p_values_sal] = bootstrapRidge(X_mat, abs(behav_ratingsP), 1000, 0.1);

    Pmat_val =X_mat(:,p_values_val<0.05);
    Pmat_sal = X_mat(:,p_values_sal<0.05);
    P_val = corrcoef(Pmat_val');
    P_val = P_val(utl_mask);
    P_sal = corrcoef(Pmat_sal');
    P_sal = P_sal(utl_mask);

    behav_ratings_ = squeeze(h2(s,:,:));
    bv = squeeze(h2_sal(s,:,:));

    M_val = corrcoef(behav_ratings_');
    M_val = M_val(utl_mask);
    M_sal = corrcoef(bv');
    M_sal = M_sal(utl_mask);

    matr(s,1) = fastcorr(M_val,P_val);
    matr(s,2) = fastcorr(M_sal,P_sal);
end
bar(matr)
yline(r2t(0.05,12720))

%% DNN performance
DNNloc = 'C:\Work\ARC\ARC_DNN\sp_val-12-3-7';
filename = fullfile(DNNloc,'variables.json'); % replace with the path to your file
fileText = fileread(filename);
dataStruct = jsondecode(fileText);

% rpred = dataStruct.r_coef;
% % Yuck
% rpred2(:,1)=rpred(1,:,1);
% rpred2(:,2)=rpred(1,:,2);
% rpred2(:,3)=rpred(2,:,1);
% rpred2(:,4)=rpred(2,:,2);
% rpred2(:,5)=rpred(3,:,1);
% rpred2(:,6)=rpred(3,:,2);
% boxplot(rpred2)


rpred = dataStruct.p_mat;
% Yuck
rpred2(:,1)=rpred(1,2,:);
rpred2(:,2)=rpred(1,3,:);
rpred2(:,3)=rpred(2,2,:);
rpred2(:,4)=rpred(2,3,:);
rpred2(:,5)=rpred(3,2,:);
rpred2(:,6)=rpred(3,3,:);
boxplot(rpred2)