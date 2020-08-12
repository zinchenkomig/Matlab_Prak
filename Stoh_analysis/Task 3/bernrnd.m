function x = bernrnd(p,varargin)
if (p>1) || (p<0)
    error('p = %4.2f is out of bounds [0, 1]',p)
end
x = rand(varargin{:}) < p;
end