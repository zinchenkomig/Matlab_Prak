general_sample_size = 1e7;
min_sample_size = 1e3;
mu = 15;
sigma = 3;
alpha = 0.9;
X = random('Normal', mu, sigma, [1, general_sample_size]);
nVec = linspace(min_sample_size, general_sample_size, 25);
sample_mean = zeros(size(nVec));
sample_variance = zeros(size(nVec));
left_mean = zeros(size(nVec));
right_mean = zeros(size(nVec));
left_var = zeros(size(nVec));
right_var = zeros(size(nVec));

for i = 1:length(nVec)
    sample_mean(i) = mean(X(1:nVec(i)));
    sample_variance(i) = var(X(1:nVec(i)));
    
    t = tinv((1 - alpha)/2, nVec(i) - 1);
    left_mean(i) = sample_mean(i) - (sample_variance(i) / nVec(i)) ^ 0.5 * t;
    right_mean(i) = sample_mean(i) + (sample_variance(i) / nVec(i)) ^ 0.5 * t;
    
    left_chi = chi2inv((1+alpha)/2, nVec(i) - 1);
    right_chi = chi2inv((1 - alpha)/2, nVec(i) - 1);
    
    left_var(i) = (nVec(i) - 1) * (sample_variance(i))/ (left_chi);
    right_var(i) = (nVec(i) - 1) * (sample_variance(i)) / (right_chi);
end

errorbar(nVec, sample_mean, sample_mean-left_mean, right_mean - sample_mean);
set(gca,'Xscale','log');
xlabel('sample size');
ylabel('\mu');
figure;

errorbar(nVec, sample_variance, sample_variance-left_var, right_var - sample_variance);
xlabel('sample size');
ylabel('\sigma^2');
set(gca,'Xscale','log');


