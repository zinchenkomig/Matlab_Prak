function cauchy = cauchyRnd(x0, gamma,varargin)
    x = rand(varargin{:});
    cauchy = x0 + gamma * tan(pi * (x - 0.5));
end