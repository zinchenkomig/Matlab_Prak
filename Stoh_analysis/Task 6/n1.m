%% Monte-Carlo
n = 1e5;
fSum = 0;
experiments_amount = 1;
epsilon = 0;
times_monte = zeros(1, experiments_amount);
disp('Monte-Carlo');
for i = 1:experiments_amount
    tstart = tic();
    XMatr = 0.5 ^ 0.5 * randn(n, 10);
    sqProd = prod(XMatr, 2) .^ 2;
    integralF = pi ^ 5 * exp(- 0.5 ^ 7 ./ sqProd) ./ sqProd;
    sigma = var(integralF) ^ (0.5);
    epsilon =  epsilon + 3 * sigma / (n ^ 0.5);
    fSum = sum(integralF) + fSum;
    times_monte(i) = toc(tstart);
end
integralMonte = fSum / (n * experiments_amount);
t = mean(times_monte);
epsilon = epsilon / experiments_amount;
fprintf('I = %4.2f \n', integralMonte);
fprintf('t = %4.4f \n', t);
fprintf('error = %4.2f \n', epsilon);


%% Rectangles
n = 5;
dimensions = 10;
xGrid = linspace(0,1,n + 2);
xGrid = xGrid(2:n+1);
N = xGrid(1);
i = zeros(10, n);
integralRect = 0;
disp('Rectangles');
disp('Please wait...');
tstart = tic();
for i1 = 1:n
    for i2 = 1:n
        for i3 = 1:n
            for i4 = 1:n
                for i5 = 1:n
                    for i6 = 1:n
                        for i7 = 1:n
                            for i8 = 1:n
                                for i9 = 1:n
                                    for i10 = 1:n
                                        z = [xGrid(i1), xGrid(i2), xGrid(i3), xGrid(i4), xGrid(i5), xGrid(i6), xGrid(i7), xGrid(i8), xGrid(i9), xGrid(i10)];
                                        tanSq = tan(pi/2 * z) .^ 2;
                                        integralRect = integralRect + (pi*N)^10 * ...
                                            exp(-(sum(tanSq) + 1./(2^7 * prod(tanSq)))) /...
                                            prod(sin(pi/2*z)) ^ 2;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
tRect = toc(tstart);

fprintf('I = %4.2f \n', integralRect);
fprintf('t = %4.4f \n', tRect);
