f = @(x)((x^3) - (4 * x^2) + (5 * x) - (2));
f1 = @(x)((3 * x^2) - (8 * x) + (5));
itmax = 10^(10);
table = zeros(10, 3);
    
for j = 1 : 10
    x0 = 0;
    tolx = 10^(-j);
    table(j,1) = tolx;
   
    i = 0;
    x = x0;
    vai = 1;
   while ((i < itmax) && vai)
        i = i + 1;
        x0 = x;
        fx = feval(f, x0);
        f1x = feval(f1, x0);
        table(j ,2) = table (j, 2) + 2;
        x1 = x0 - (fx / f1x);
        fx = feval(f, x1);
        f1x = feval(f1, x1);
        table(j ,2) = table (j, 2) + 2;
        x = x1 - (fx / f1x);
        if (x - (2 * x1) + x0) ~= 0
            x = ((x * x0) - (x1^2)) / (x - (2 * x1) + x0);
            table(j, 3) = table(j, 3) + 1;
        else
            table(j, 3) = i;
            i = itmax;
        end
        vai = abs(x - x0) > tolx;
    end
    if vai
        disp('il metodo non converge');
    end
end