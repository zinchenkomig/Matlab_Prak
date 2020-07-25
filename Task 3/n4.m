n = 1e4;
alpha = 0.05;
nTests = 100;
accepted_e = 0;
accepted_v = 0;
for i = 1:nTests
    ksi = standard_pairs(n);
    
    t = mean(ksi) / ((var(ksi) / n)^(0.5));
    critSt = tinv(1 - alpha/2, n-1);
    if abs(t)<critSt
        accepted_e = accepted_e + 1;
    end
    
    normVec = randn(size(ksi));
    F = var(ksi) / var(normVec);
    fleft = finv(alpha/2, n-1, n-1);
    fright = finv(1 - alpha/2, n-1, n-1);
    if F>fleft && F<fright
        accepted_v = accepted_v + 1;
    end
end

fprintf('H0 (Expectance) Accepted: %4.2f %% \n', 100 * accepted_e/nTests);
fprintf('H0 (Variance) Accepted: %4.2f %% \n', 100 * accepted_v/nTests);


histogram(ksi, 'Normalization', 'pdf');
hold on;
xVec = linspace(-3,3);
plot(xVec, normpdf(xVec), 'Color', [0.9, .2, .2], 'LineWidth', 1.3);
legend('Generated', 'Real PDF');
xlabel('random variable');
ylabel('probability');







