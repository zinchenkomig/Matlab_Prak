function [x,F] = cantrnd(n,eps)
    if nargin < 1
        n = 1; % объем выборки 
    end
    if nargin < 2
        eps = 1e-10; % точность расчета
    end
    m = round(-log(eps)/log(3));
    bern = bernrnd(0.5,n,m);
    deg = -(1:m)';
    x = 2 * bern * 3.^deg;
    F = bern * 2.^deg;
end
