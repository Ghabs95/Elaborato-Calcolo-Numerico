f = @(x)((x^3) - (4 * x^2) + (5 * x) - (2));
f1 = @(x)((3 * x^2) - (8 * x) + (5));
itmax = 10^(10);
table = zeros(10, 3);
    
for j = 1 : 10
    x0 = 0;
    mol = 2;
    tolx = 10^(-j);
    table(j,1) = tolx;
        
    fx = feval(f, x0);
    f1x = feval(f1, x0);
    table(j ,2) = table (j, 2) + 2;
    x = x0 - (fx / f1x);
    i = 0;
    while((i < itmax) && (abs(x - x0) > tolx))
        i = i + 1;
        x0 = x;
        fx = feval(f, x0);
        f1x = feval(f1, x0);
        table(j ,2) = table (j, 2) + 2;
        x = x0 - mol * (fx / f1x);
    end
    table(j ,3) = i + 1;
    if (abs(x - x0) > tolx)
        disp('il metodo non converge');
    end
end