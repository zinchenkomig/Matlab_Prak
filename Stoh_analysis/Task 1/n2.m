%Calculations
m = 3;
geomVec = geometric(0.3, 10000, 1);

%Plotting
figure;
histogram(geomVec, 'Normalization', 'probability');
hold on;
plot((0:30), geopdf(0:30, 0.3), 'Color', [0.9, 0.2, 0.2], 'LineWidth', 1.3);
legend('Computed', 'Real PDF');
xlabel('random variable');
ylabel('probability');

figure;
histogram(geomVec, 'Normalization', 'probability');
hold on;
histogram(geomVec(geomVec>m) - m - 1, 'Normalization', 'probability');
legend('P(X > m)', 'P(X > m + n| X > n)');
xlabel('random variable');
ylabel('probability');