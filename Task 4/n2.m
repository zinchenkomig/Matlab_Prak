n = 1e4;
normalVec = standard_phoneNewMan(n);
% histogram(normalVec, 'Normalization', 'pdf');
mu = 3;
sigma = 1;
normalVec1 = mu + sigma*normalVec;
mu = 0;
sigma = 4;
normalVec2 = mu + sigma*normalVec;
normalVec = [normalVec',normalVec1', normalVec2'];
normplot(normalVec);
xlabel('random variable');
legend('\mu = 0; \sigma = 1', '\mu = 3; \sigma = 1', '\mu = 0; \sigma = 4');
