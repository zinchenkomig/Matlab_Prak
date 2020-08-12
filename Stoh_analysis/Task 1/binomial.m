function x = binomial(p, n, varargin)
x = sum(bernrnd(p, n, varargin{:}), 1);
x = reshape(x, varargin{:});
end