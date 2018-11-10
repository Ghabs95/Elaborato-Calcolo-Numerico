 x0 = 5;
 x1 = 3;
 tolx = 10^-10;
 itmax = 10^10;
 a = 5;

f = @(x)((x^2) - a);

table = zeros;

fx0 = feval(f, x0);
fx1 = feval(f, x1);
x = ((fx1 * x0) - (fx0 * x1))/(fx1-fx0);
i = 0;
table(i + 1, 1) = x;
table(i + 1, 2) = abs((a^(1/2))-x)/(a^(1/2));
while((i < itmax) && (abs(x - x1) > tolx))
    i = i + 1;
    x0 = x1;
    x1 = x;
    fx0 = feval(f, x0);
    fx1 = feval(f, x1);
    x = ((fx1 * x0) - (fx0 * x1))/(fx1-fx0);
    table(i + 1, 1) = x;
    table(i + 1, 2) = abs((a^(1/2))-x)/(a^(1/2));
end
if (abs(x - x1) > tolx)
    disp('il metodo non converge');
end