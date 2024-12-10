% Set seed for reproducibility
rng(123);

% Number of samples and instances
N = 160;
instances = 10;

% Generate mean vectors for x and y
x_mean = randn(N, 1); 
y_mean = 0.1 * x_mean + randn(N, 1); % y is weakly correlated with x

% Generate x and y samples based on their means
x_samples = zeros(N, instances);
y_samples = zeros(N, instances);
for i = 1:instances
    x_samples(:, i) = x_mean + 0.5 * randn(N, 1);
    y_samples(:, i) = y_mean + 0.5 * randn(N, 1);
end

% Compute correlation between x_mean and y_mean
corr_mean = corr(x_mean, y_mean);

% Compute average correlation between individual instances
correlations = zeros(instances, instances);
for i = 1:instances
    for j = i+1:instances
        correlations(i, j) = corr(x_samples(:, i), y_samples(:, j));
    end
end
avg_corr_instances = mean(correlations(correlations~=0));

disp(['Correlation between x_mean and y_mean: ', num2str(corr_mean)]);
disp(['Average correlation between individual instances of x and y: ', num2str(avg_corr_instances)]);
