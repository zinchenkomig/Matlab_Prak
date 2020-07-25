n = 1e3;
a = 1;
b = 0;
samples = 1e4;
x = cauchyRnd(b, a,samples, n);
sn = cumsum(x, 2);
nVec = 1:n;
xVec = sn(1,:) ./ nVec;

plot(nVec, xVec);
xlabel('n');
ylabel('S_n/n');
legend('S_n/n');

figure;
histogram(sn(:,end) / n, 'Normalization', 'pdf', 'BinWidth',0.2);
xVec = linspace(-5, 5);
xlim(5 .* [-1, 1])
y = CauchyPdf(a,b,xVec);
hold on;
plot(xVec, y, 'Color', [0.9, .2, .2], 'LineWidth', 1.3);
xlabel('random variable');
ylabel('probability');
legend('S_n/n', 'Cauchy PDF');


function y = CauchyPdf(a,b,x)
y = 1 ./ (pi * a * (1 + ((x - b)/ a).^2));
end
