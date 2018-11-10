function If = trapcomp(n, a, b, fun)
% If = trapcomp(n, a, b, fun)
% Calcola l'integrale della funzione, nell'intervallo prescelto,
% usando la formula dei trapezi composita.
	if n <= 0, error('Numero intervalli non valido.'), end
    if a >= b, error('Intervallo di integrazione non valido.'), end
	x = linspace(a, b, n+1);
	f = feval(fun, x);
	If = ((b-a)/n)*(f(1)/2 + sum(f(2:n)) + f(n+1)/2);
end

function If = simpcomp(n, a, b, fun)
% If = simpcomp(n, a, b, fun)
% Calcola l'integrale della funzione, nell'intervallo prescelto,
% usando la formula di Simpson composita.
    if mod(n,2) ~= 0
		error('Il numero di intervalli deve essere pari.')
    end
    if n <= 0
        error('Numero di intervalli non valido.')
    end
    if a >= b, error('Intervallo di integrazione non valido.'), end
	x = linspace(a, b, n+1);
	f = feval(fun, x);
	If = f(1) + f(n+1);
	If = (If + 4*sum(f(2:2:n)) + 2*sum(f(3:2:n-1)))*(b-a)/(3*n);
end

function If = trapad(a, b, fun, tol, fa, fb)
% If = trapad(a, b, fun, tol)
% Calcola ricorsivamente l'integrale della funzione, nell'intervallo prescelto, 
% usando la formula dei trapezi adattiva.
	if a >= b, error('Intervallo di integrazione non valido.'), end
	if nargin <= 4
		fa = feval(fun, a);
		fb = feval(fun, b);
	end
	h = b - a;
	x1 = (a+b)/2;
	f1 = feval(fun, x1);
	I1 = (h/2)*(fa + fb);
	If = (I1 + h*f1)/2;
	err = abs(If - I1)/3;
    if err > tol
		If = trapad(a, x1, fun, tol/2, fa, f1) ...
           + trapad(x1, b, fun, tol/2, f1, fb);
    end
end

function If = simpad(a, b, fun, tol, fa, f1, fb)
% If = simpad(a, b, fun, tol)
% Calcola ricorsivamente l'integrale della funzione, nell'intervallo prescelto, 
% usando la formula di Simpson adattiva.
	if a >= b, error('Intervallo di integrazione non valido.'), end
	x1 = (a+b)/2;
    if nargin <= 4
		fa = feval(fun, a);
		fb = feval(fun, b);
		f1 = feval(fun, x1);
    end
    h = (b - a)/6;
    I1 = h*(fa + 4*f1 + fb);
	f2 = feval(fun, (a+x1)/2);
	f3 = feval(fun, (x1+b)/2);
	If = .5*h*(fa + 4*f2 + 2*f1 + 4*f3 + fb);
	err = abs(If - I1)/15;
    if err > tol
		If = simpad(a, x1, fun, tol/2, fa, f2, f1)...
            + simpad(x1, b, fun, tol/2, f1, f3, fb);
    end
end
