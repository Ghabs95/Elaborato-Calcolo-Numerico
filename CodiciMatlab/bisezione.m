function x = bisezione(f, a, b, tolx)
% x = bisezione(f, a, b, tolx)
% Questa funzione implementa il metodo di bisezione sulla funzione f
% nell'intervallo a, b con tolleranza tolx.
% f = funzione di cui trovale lo zero
% a = primo estremo dell'intervallo in cui cercare lo zero
% b = secondo estremo dell'intervallo in cui cercare lo zero
% tolx = errore con cui calcolare lo zero della funzione
    fa = feval(f, a);
    fb = feval(f, b);
    x = (a + b) / 2;
    fx = feval(f, x);
    imax = ceil(log2(b - a) - log2(tolx));
    for i = 2 : imax
        f1x = abs((fb - fa) / (b - a));
        if abs(fx) <= tolx * f1x
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
    end