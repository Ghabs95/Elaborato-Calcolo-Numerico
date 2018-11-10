function x = root_newton(x0, tolx, itmax, a)
% x = root_newton(x0, tolx, itmax, a)
% Questo metodo approssima la radice di a partendo da x0 usando il metodo
% di newton con tolleranza tolx e numero massimo di iterazioni itmax.
% x0 = punto di innesto del metodo
% tolx = tolleranza nell'approssimazione dello zero di f
% itmax = numero massimo di iterazioni 
% a = numero di cui approssimare la radice

   
    f = @(x)((x^2) - a);
    f1 = @(x)(2 * x);
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