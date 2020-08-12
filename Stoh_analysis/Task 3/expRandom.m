function x = expRandom(lambda, varargin)
x = - log(rand(varargin{:})) ./ lambda;
end