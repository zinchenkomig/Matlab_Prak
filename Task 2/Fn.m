function y = Fn(xVec, xSample)
n = size(xSample, 1);
if size(xVec, 1) == 1
    xVec = xVec';
end
y = (1/n) * sum(xSample <= xVec')';
end