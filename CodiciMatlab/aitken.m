function x = aitken(f, f1, x0, tolx, itmax)
% x = aitken(f, f1, x0, tolx ,itmax)
% La funzione implementa il metodo dell'accelerazione di aitken sulla
% funzione f dal punto x0 con tolleranza tolx ed in un massimo di itmax
% iterazioni.
% f = funzione della quale trovare lo zero
% f1 = derivata prima di f
% x0 = punto di innesto del metodo
% tolx = tolleranza nell'approssimazione dello zero di f
% itmax = numero massimo di iterazioni 

    i = 0;
    x = x0;
    vai = 1;
    while ((i < itmax) && vai)
        i = i + 1;
        x0 = x;
        fx = feval(f, x0);
        f1x = feval(f1, x0);
        x1 = x0 - (fx / f1x);
        fx = feval(f, x1);
        f1x = feval(f1, x1);
        x = x1 - (fx / f1x);
        if (x - (2 * x1) + x0) ~= 0
            x = ((x * x0) - (x1^2)) / (x - (2 * x1) + x0);
        else 
            i = itmax;
        end
        vai = abs(x - x0) > tolx;
    end
    if vai
        disp('il metodo non converge');
    end