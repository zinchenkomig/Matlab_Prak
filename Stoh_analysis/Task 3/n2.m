lambda = 10;
n = 1e5;
x = PoisRnd(lambda, n, 1);
histogram(x, 'Normalization', 'probability');
hold on;
plot((0:25), poisspdf((0:25), lambda), 'Color', [0.9, .2, .2], 'LineWidth', 1.3);
legend('Generated', 'Real PDF');
xlabel('random variable');
ylabel('probability');


function kMatr = PoisRnd(lambda, varargin)
sumsMatr = zeros(varargin{:});
kMatr = zeros(varargin{:});
while any(sumsMatr < 1)
    x = expRandom(lambda, varargin{:});
    sumsMatr = sumsMatr + x;
    kMatr(sumsMatr < 1) = kMatr(sumsMatr < 1) + 1; 
end
end