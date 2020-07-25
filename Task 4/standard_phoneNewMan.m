function normalVec = standard_phoneNewMan(n)
left = n;
first = 1;
normalVec = zeros(1,n);
while left ~= 0
    xVec = cauchyRnd(0,1,1,left);
    pBern = (exp(0.5))/2 * exp(-xVec.^2 / 2) .* (xVec.^2 + 1);
    vVec = bernrnd(pBern, 1, left);
    used = sum(vVec);
    left = left - used;
    normalVec(first:first+used-1) = xVec(vVec == 1);
    first = first + used;
end
end