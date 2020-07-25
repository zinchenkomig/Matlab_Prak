n = 1e5;
p = 0.99;
tstart = tic();
phi = rand(1, n) * 2 * pi;
r = rand(1, n);

x1 = r.^0.5 .* cos(phi);
x2 = r.^0.5 .* sin(phi);

y = min(f(x1,x2));
t = toc(tstart);

epsilon = 52 * (p / n)^0.5;

disp('Random Search');
fprintf('min f = %4.4f \n', y);
fprintf('eps = %4.2f \n', epsilon);
fprintf('t = %4.4f \n', t);

function y = f(x1,x2)
y = x1 .^ 3 .* sin(1 ./ x1) + 10 * x1 .* x2 .^ 4 .* cos(1 ./ x2);
end
