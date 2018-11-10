f = @(x)((x^3) - (4 * x^2) + (5 * x) - (2));
f1 = @(x)((3 * x^2) - (8 * x) + (5));
itmax = 10^(10);
table = zeros(20, 3);
    
for j = 1 : 20
    x1 = 3;
    tolx = 10^(-j);
    table(j,1) = tolx;
        
    fx = feval(f, x1);
    f1x = feval(f1, x1);
    table(j ,2) = table (j, 2) + 2;
    x = x1 - (fx / f1x);
    i = 0;
    while((i < itmax) && (abs(x - x1) > tolx))
        i = i + 1;
        x0 = x1;
        x1 = x;
        fx0 = feval(f, x0);
        fx1 = feval(f, x1);
        table(j ,2) = table (j, 2) + 2;
      x = ((fx1 * x0) - (fx0 * x1))/(fx1-fx0);
    end
    table(j ,3) = i + 1;
    if (abs(x - x1) > tolx)
        disp('il metodo non converge');
    end
end