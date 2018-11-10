function x = newton(f, f1, x0, tolx, itmax)
% x = newton(f, f1, x0, tolx, itmax)
% Questa funzione applica il metodo di newton per trovare uno zero della
% funzione f di derivata prima f1 a partire dal punto x0.
% f = funzione della quale trovare lo zero
% f1 = derivata prima di f
% x0 = punto di innesto del metodo
% tolx = tolleranza nell'approssimazione dello zero di f
% itmax = numero massimo di iterazioni 

    fx = feval(f, x0);
    f1x = feval(f1, x0);
    x = x0 - (fx / f1x);
    i = 0;
    while((i < itmax) && (abs(x - x0) > tolx))
        i = i + 1;
        x0 = x;
        fx = feval(f, x0);
        f1x = feval(f1, x0);
        x = x0 - (fx / f1x);
    end
    if (abs(x - x0) > tolx)
        disp('il metodo non converge');
    end