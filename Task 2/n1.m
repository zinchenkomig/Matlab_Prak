n = 1e4;
experiments = 20;
alpha = 0.05;
positive_count = 0;

for i = 1:experiments
    fprintf('%d of %d \n',i,experiments);
    [xSample,F] = cantrnd(n);
    if (kolmCrit(xSample, F, alpha))
        positive_count = positive_count + 1;
    end
end
fprintf('H0 Accepted: %4.2f %% \n', 100 * positive_count/experiments);

%Plotting
F = sort(F);
xSample = sort(xSample);
plot(xSample,F);
xlabel('x');
ylabel('F');

function y = kolmCrit(xSample, F, alpha)

Dn = max(abs(Fn(xSample,xSample) - F));
n = size(xSample,1);
K_alpha = fsolve(@(x) kolmcdf(x) - (1 - alpha), 1);

if(n^0.5 * Dn < K_alpha)
    y = 1;
else
    y = 0;
end

end
