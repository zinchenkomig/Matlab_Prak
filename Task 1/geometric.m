function x = geometric(p,varargin)
x = zeros(varargin{:});
SuccessMatr = bernrnd(p,varargin{:});
while any(~SuccessMatr)
    x(~SuccessMatr) = x(~SuccessMatr) + 1;
    SuccessMatr = SuccessMatr | bernrnd(p,varargin{:});
end
end