function x = root_secanti( x0, x1, tolx, itmax, a)
% x = root_secanti(x0, x1, tolx, itmax, a)
% Questo metodo approssima la radice di a partendo da x0 e da x1 usando il 
% metodo delle secanti con tolleranza tolx e numero massimo di iterazioni 
% itmax.
% x0 = primo punto di innesto del metodo
% x1 = secondo punto di innesto del metodo
% tolx = tolleranza nell'approssimazione dello zero di f
% itmax = numero massimo di iterazioni 
% a = numero di cui approssimare la radice

 
    f = @(x)((x^2) - a);
    fx0 = feval(f, x0);
    fx1 = feval(f, x1);
    x = ((fx1 * x0) - (fx0 * x1))/(fx1-fx0);
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