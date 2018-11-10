f = @(x) -2*(x^(-3))*cos(x^(-2));
If = simpad(1/2,100,f,10^(-4))

function If = simpad(a, b, fun, tol, fa, f1, fb)
% If = simpad(a, b, fun, tol)
% Calcola ricorsivamente l'integrale della funzione, nell'intervallo prescelto, 
% usando la formula di Simpson adattiva.
% Input:
 % 	a: estremo sinistro
 % 	b: estremo destro
 % 	fun: funzion integranda
 %	tol: tolleranza
% Output:
 % 	If: approssimazione dell'integrale definito della funzione
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