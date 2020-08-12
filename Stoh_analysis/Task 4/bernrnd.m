function x = bernrnd(p,varargin)
if any((p>1) | (p<0))
    error('p is out of bounds [0, 1]')
end
x = rand(varargin{:}) < p;
end
