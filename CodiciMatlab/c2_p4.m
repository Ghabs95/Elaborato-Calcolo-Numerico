x0 = 5;
a = 5;
tolx = 10^-10;
itmax = 10^10;

f = @(x)((x^2) - a);
f1 = @(x)(2 * x);

table = zeros;

fx = feval(f, x0);
f1x = feval(f1, x0);
x = x0 - (fx / f1x);
i = 0;
table(i + 1, 1) = x;
table(i + 1, 2) = abs((a^(1/2))-x)/(a^(1/2));
while((i < itmax) && (abs(x - x0) > tolx))
    i = i + 1;
    x0 = x;
    fx = feval(f, x0);
    f1x = feval(f1, x0);
    x = x0 - (fx / f1x);
    table(i + 1, 1) = x;
    table(i + 1, 2) = abs((a^(1/2))-x)/(a^(1/2));
end
if (abs(x - x0) > tolx)
    disp('il metodo non converge');
end