n0 = 1000;
n1 = 1000000;
nVec = (n0:100:n1)';
[xSample, F] = cantrnd(n1);
sumsVec = cumsum(xSample);
empE = sumsVec(nVec) ./ nVec;
sqSumsVec = cumsum(xSample .^ 2);
empEsq = sqSumsVec(nVec) ./ nVec; 
empVar = empEsq - empE .^ 2;

figure;
plot(nVec, empE, 'DisplayName', 'Emperical Expectance');
set(gca,'Xscale','log');
xlabel('n');
ylabel('Expectance');
hold on;
plot([n0, n1], [0.5, 0.5], 'DisplayName', 'Theoretical Expectance');
legend;

figure;
plot(nVec, empVar, 'DisplayName', 'Emperical Variance');
xlabel('n');
ylabel('Variance');
hold on;
plot([n0 n1], [0.125 0.125], 'DisplayName', 'Theoretical Variance');
legend;
set(gca,'Xscale','log');
