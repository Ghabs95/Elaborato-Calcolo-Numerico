table = zeros(10,3);   
for j = 1 : 10
    a = 0;
    b = 3;
    tolx = (10^(-j));
    table(j,1) = tolx;
    table(j, 2) = 0;
    fa = feval(f, a);
    fb = feval(f, b);
    table(j, 2) = table(j, 2) + 2;
    x = (a + b) / 2;
    fx = feval(f, x);
    table(j, 2) = table(j, 2)+1;
    imax = ceil(log2(b - a) - log2(tolx));
    for i = 2 : imax
        fix = abs((fb - fa) / (b - a));
        if abs(fx) <= tolx * fix
            table(j, 3) = i;
            break
        elseif fa * fx < 0
            b = x;
            fb = fx;
        else
            a = x;
            fa = fx;
        end
        x = (a + b) / 2;
        fx = feval(f, x);
        table(j, 2) = table(j, 2) + 1;
        table(j, 3) = i;
    end
end