function Fk = kolmcdf(x)
n = 1e5;
k = 1:n;
Fk = 1 - 2 * sum((-1).^(k-1) .* exp(-2*k.^2 .* x.^2), 2);
if (x<= 0)
    Fk = 0;
end
end