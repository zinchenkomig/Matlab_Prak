x1 = -210;
x2 = 154;
g_current = g(x1,x2);
q = 1;
i = 1;
while(q > 1e-7)
    x1New = tan(pi / 2 * (rand() * 2 - 1));
    x2New = tan(pi / 2 * (rand() * 2 - 1));
    g_new = g(x1New, x2New);
    if (g_new < g_current) || (rand() < exp( - (g_new - g_current) / q))
        g_current = g_new;
        x1 = x1New;
        x2 = x2New;
    end 
    i = i + 1;
    q = 1 / i;
end
disp(g_current);
disp(x1);
disp(x2);
disp(i);

function y = g(x1, x2)
y = (x1 - 1).^2 + 100 * (x2 - x1 .^ 2) .^ 2;
end