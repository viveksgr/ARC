function pilat = nii_extract_list(dirs,anat_cell,anat_names)


nS = 3;
nf = 2;
nanat = length(anat_names);
valmat = zeros(nanat,nf,nS);
pdir = fullfile(dirs,{'sal','val'});
prefixes = {'P_','P_'};


for jj = 1:nf
    pdirgroup = pdir{jj};
    for ss = 1:nS       
        for aa = 1:nanat
            flist = dir(fullfile(pdir{jj},sprintf('%sARC%02d_*_RSAx1%s.nii',prefixes{jj},ss,anat_names{aa})));
            assert(length(flist)==1)
            fmat = spm_read_vols(spm_vol(fullfile(pdir{jj},flist.name)));
            pvox = fmat(logical(squeeze(anat_cell{ss}(:,:,:,aa))));
            % valmat(aa,jj,ss)= tanh(mean(atanh(pvox)));
            valmat(aa,jj,ss)= sum(pvox>0.01)/length(pvox);
        end
    end
end

pilat  = makepilat(valmat);
xticks(1:length(anat_names))
xtickangle(90)
ylabel('mean (r)')
xticklabels(anat_names)
yline(r2t(0.05,12720))
legend({'sal','val','S1','S2','S3','p=0.05'})
savefig('datafig')
print('datafig','-dpng')

% legend(fname)
end

function pilat = makepilat(dataMat)
% Assuming data matrix is named dataMat with dimensions AXGXS
[A, G, S] = size(dataMat);

% Calculate means and standard errors over S
meanVals = mean(dataMat, 3);
stdErrs = std(dataMat, 0, 3) / sqrt(S);  % Change to std(dataMat, 0, 3) for standard deviation

% Create the grouped bar plot
pilat = figure;
barHandle = bar(meanVals);
hold on;

% Add error bars
groupwidth = min(0.8, G/(G+1.5));
for i = 1:G
    % Calculate center of each bar group
    x = (1:A) - groupwidth/2 + (2*i-1) * groupwidth / (2*G);
    errorbar(x, meanVals(:, i), stdErrs(:, i), 'k', 'linestyle', 'none','HandleVisibility','off');
end

% Add individual data points for each S
% colors = [1 0 0; 0 1 0; 0 0 1]
colors = lines(S+2);  % Generates a matrix of RGB values for S different colors
for s = 1:S
    for a = 1:A
        y = squeeze(dataMat(a, :, s));
        x = a - groupwidth/2 + (2*(1:G)-1) * groupwidth / (2*G);
        if a==1
            plot(x, y, '-o', 'Color', colors(s+2, :), 'MarkerFaceColor', colors(s+2, :));
        else
            plot(x, y, '-o', 'Color', colors(s+2, :), 'MarkerFaceColor', colors(s+2, :),'HandleVisibility','off');
        end
    end
end

end