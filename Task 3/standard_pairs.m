function ksi = standard_pairs(n)
n = floor(n/2);
w = expRandom(0.5, 1, n);
phi = 2 * pi * rand(1, n);
ksi = [w.^(1/2) .* cos(phi), w.^(1/2) .* sin(phi)];
end