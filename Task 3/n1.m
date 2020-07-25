n = 1e5;
lambda1 = 0.2;
lambda2 = 0.4;
lambda3 = 0.5;
lambda = lambda1;
m = size(lambda, 1);
t = 2;
alpha = 0.005;

x1 = expRandom(lambda1, n, 1);

x2 = expRandom(lambda2, n, 1);
x3 = expRandom(lambda3, n, 1);

xMin = min([x1, x2, x3]');
xTheor = expRandom(lambda1 + lambda2 + lambda3, n, 1);

figure;
histogram(x1, 'Normalization', 'probability','BinWidth',1);
hold on;
plot((0:30), exppdf(0:30, 1/lambda1), 'Color', [0.9, .2, .2], 'LineWidth', 1.3);
legend('Generated', 'Real PDF')
xlabel('random variable');
ylabel('probability');

figure;
histogram(x1, 'Normalization', 'probability', 'DisplayName', 'P(Y>m)','BinWidth',0.7);
hold on;
histogram(x1(x1 > t) - t, 'Normalization', 'probability','DisplayName', 'P(Y>n+m|Y>n)','BinWidth',0.7);
xlim([0 invF(1 - alpha, lambda1) ]);
xlabel('random variable');
ylabel('probability');
legend;

figure;
histogram(xMin, 'Normalization', 'probability', 'DisplayName', 'Emperical','BinWidth',0.2);
hold on;
histogram(xTheor, 'Normalization', 'probability', 'DisplayName', 'Theoretical','BinWidth',0.2);
xlabel('random variable');
ylabel('probability');
legend;





