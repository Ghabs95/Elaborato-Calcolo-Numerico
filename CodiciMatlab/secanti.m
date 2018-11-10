function x = secanti(f, f1, x1, tolx, itmax)
% x = secanti(f, f1, x1, tolx, itmax)
% La funzione applica il metodo delle secanti alla funzione f con punto di
% innesto x1, tolleranza tolx enumero massimo di iterazioni itmax
% f = funzione della quale trovare lo zero
% f1 = derivata prima di f
% x1 = punto di innesto del metodo
% tolx = tolleranza nell'approssimazione dello zero di f
% itmax = numero massimo di iterazioni 

    fx = feval(f, x1);
    f1x = feval(f1, x1);
    x = x1 - (fx / f1x);
    i = 0;
    while((i < itmax) && (abs(x - x1) > tolx))
        i = i + 1;
        x0 = x1;
        x1 = x;
        fx0 = feval(f, x0);
        fx1 = feval(f, x1);
        x = ((fx1 * x0) - (fx0 * x1))/(fx1-fx0);
    end
    if (abs(x - x1) > tolx)
        disp('il metodo non converge');
    end