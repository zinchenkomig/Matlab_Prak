%Calculations
binomialVec = binomial(0.3, 100, 10000, 1);

%Plotting
figure;
histogram(binomialVec, 'Normalization', 'probability');
hold on;
plot((10:50), binopdf(10:50, 100, 0.3), 'Color', [0.9, .2, .2], 'LineWidth', 1.3);
xlabel('random variable');
ylabel('probability');
legend('Computed', 'Real PDF');
