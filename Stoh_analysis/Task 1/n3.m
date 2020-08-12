%Calculations
n = 100000;
x = bernrnd(0.5, n, 1);
c = x + (x - 1);
nVec = (1:n)';
xVec = cumsum(c) ./ (n .^ 0.5);

%Plotting
plot(xVec);
xlabel('n');
ylabel('Y(n)');