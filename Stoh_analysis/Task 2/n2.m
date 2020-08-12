%Calculations
n = 1e4;
experiments = 20;
alpha = 0.05;
symmetry_count = 0;
onethird_count = 0;
for i = 1:experiments
    fprintf('%d of %d \n',i,experiments);
    [xSample, F] = cantrnd(n);
    y = SmirnovCrit(xSample, 1 - xSample, alpha);
    if y <= 0
        symmetry_count = symmetry_count + 1;
        %disp('H0 (x, 1 - x): Accepted');
    else
        %disp('H0 (x, 1 - x): Declined');
    end
    y = SmirnovCrit(xSample/3, xSample(xSample < (1/3)), alpha);
    if y <= 0
        onethird_count = onethird_count + 1;
        %disp('H0 (x/3, x in [0, 1/3]): Accepted');
    else
        %disp('H0 (x/3, x in [0, 1/3]): Declined');
    end
end

fprintf('H0 (x, 1 - x) Accepted: %4.2f %% \n', 100 * symmetry_count/experiments);
fprintf('H0 (x/3, x in [0, 1/3]) Accepted: %4.2f %% \n', 100 * onethird_count / experiments);

%Plotting
xSample = sort(xSample);
F = sort(F);

%Symmetry
figure('Name', 'Symmetry', 'NumberTitle', 'off');
plot(xSample,F);
hold on;
plot(1 - xSample, 1 - F + 3e-3,'r');
legend('x', '1 - x');
xlabel('x');
ylabel('F(x)');

%One third
figure('Name', 'One third', 'NumberTitle', 'off');
plot(xSample/3, F/2);
hold on;
plot(xSample(xSample < (1/3)), F(xSample<1/3)+ 3e-3, 'r');
legend('x/3', 'x \in [0, 1/3]');
xlabel('x');
ylabel('F(x)');

function y = SmirnovCrit(xSample, ySample, alpha)
n = max([size(xSample, 1), size(ySample,1)]);
m = min([size(xSample, 1), size(ySample,1)]);
xVec = linspace(0,1,n);
Dnm = max(abs(Fn(xVec, xSample) - Fn(xVec, ySample)));
y = Dnm - c(alpha) * ((n + m)/(n*m)) ^ 0.5;

end

