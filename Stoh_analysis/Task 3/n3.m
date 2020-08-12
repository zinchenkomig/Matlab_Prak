n = 200;
alpha = 0.05;
nTests = 100;
accepted = 0;
sampleSize = 1e4;
lambda = 4;
p = lambda / n;

for i = 1:nTests
x = binomial(p, n, sampleSize, 1);
xVec = 0:max(x)+1;

nVec = histcounts(x, xVec);
kVec = (0:max(x));
poisP = lambda .^ kVec * exp(-lambda) ./ factorial(kVec);
hi2 = sampleSize * sum((nVec / sampleSize - poisP).^2 ...
                                    ./ poisP);
Crit = chi2inv(1 - alpha, max(x));
if hi2 < Crit
    accepted = accepted + 1;
end

end

fprintf('H0 (Pois) Accepted: %4.2f %% \n', 100 * accepted/nTests);

histogram(x, 'Normalization', 'probability');
hold on;
plot((0:12), poisspdf((0:12), lambda), 'Color', [0.9, .2, .2], 'LineWidth', 1.3);
legend('Generated', 'Real PDF');
xlabel('random variable');
ylabel('probability');
