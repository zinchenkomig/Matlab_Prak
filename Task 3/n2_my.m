lambda = 10;
n = 10000;
x = PoisRnd(lambda, n, 1);
histogram(x, 'Normalization', 'probability');


function k = PoisRnd(lambda, varargin)
x = expRandom(lambda, varargin{:});
k = zeros(varargin{:});
Series = exp(-lambda);
xi = (log(1 - Series) / (-lambda));
SMatr = x > xi;
n = 0;
while any(SMatr)
    k(SMatr) = k(SMatr) + 1;
    n = n + 1;
    Series = Series + (lambda^n * exp(-lambda) / factorial(n));
    xi = (log(1 - Series) / (-lambda));
    SMatr = x>xi;
end
end