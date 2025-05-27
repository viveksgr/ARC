% Debug
figure('Position', [0.5 0.5 1280 240])
for ss = 1:3
    subplot(2,3,ss)
    hold on
    ARC_scatterDensity(t_score_mat1{ss,3}(:,1),t_score_mat2{ss,3}(:,1))
    xlabel('Old RSA')
    ylabel('New RSA')
    title(sprintf('RSA val+ OFC Sub: %02d',ss))

    subplot(2,3,ss+3)
    hold on
    ARC_scatterDensity(t_score_mat1{ss,3}(:,2),t_score_mat2{ss,3}(:,2))
    xlabel('Old RSA')
    ylabel('New RSA')
    title(sprintf('RSA val- OFC Sub: %02d',ss))
end


thr1 = 1.98;
thr2 = 1.98;
figure('Position', [0.5 0.5 1280 240])
edge = [0:0.5:20];
for ss = 1:3
    thr = tinv(0.90,size(t_score_mat1{s,3},1));
    pospop1 = and(t_score_mat1{s,3}(:,1)>thr1,t_score_mat1{s,3}(:,2)<thr);
    negpop1 = and(t_score_mat1{s,3}(:,2)>thr2,t_score_mat1{s,3}(:,1)<thr);

    % thr1 = 1.65;
    % thr2 = 1.65;
    % pospop2 = and(t_score_mat2{s,3}(:,1)>thr1,t_score_mat2{s,3}(:,2)<thr);
    % negpop2 = and(t_score_mat2{s,3}(:,2)>thr2,t_score_mat2{s,3}(:,1)<thr);

    subplot(2,3,ss)
    hold on
    histogram(t_score_mat1{s,3}(pospop1,1),edge)
    histogram(t_score_mat2{s,3}(pospop1,1),edge)
    title(sprintf('RSA val+ OFC Sub: %02d',ss))

    subplot(2,3,ss+3)
    hold on
    histogram(t_score_mat1{s,3}(negpop1,2),edge)
    histogram(t_score_mat2{s,3}(negpop1,2),edge)
    title(sprintf('RSA val- OFC Sub: %02d',ss))
end