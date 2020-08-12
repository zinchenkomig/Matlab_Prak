n = 1e5;
x0 = 0;
gamma = 1;
Cauchy = cauchyRnd(x0, gamma, 1, n);
edges = linspace(-10 * gamma, 10*gamma, 80);

histogram(Cauchy, edges, 'Normalization', 'pdf');
xlim(x0 + 9 * gamma * [-1, 1]);
legend('Cauchy distribution');
xlabel('random variable');
ylabel('probability');