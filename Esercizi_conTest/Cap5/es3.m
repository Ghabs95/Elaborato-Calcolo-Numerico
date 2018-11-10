f = @(x) -2*(x^(-3))*cos(x^(-2));
If = trapad(1/2,100,f,10^(-4))

function If = trapad(a, b, fun, tol, fa, fb)
% If = trapad(a, b, fun, tol)
% Calcola ricorsivamente l'integrale della funzione, nell'intervallo prescelto, 
% usando la formula dei trapezi adattiva.
% Input:
 % 	a: estremo sinistro
 % 	b: estremo destro
 % 	fun: funzione integranda
 %	tol: tolleranza
% Output:
 % 	If: approssimazione dell'integrale definito della funzione
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