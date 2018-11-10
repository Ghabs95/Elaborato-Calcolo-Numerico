function If = trapcomp(n, a, b, fun)
% If = trapcomp(n, a, b, fun)
% Calcola l'integrale della funzione, nell'intervallo prescelto,
% usando la formula dei trapezi composita.
% Input:
 % 	n: intero positivo, indica num intervalli in [a,b]
 % 	a: estremo sinistro
 % 	b: estremo destro
 % 	fun: funzione integranda
% Output:
 % 	If: approssimazione dell'integrale definito della funzione
	if n <= 0, error('Numero intervalli non valido.'), end
    if a >= b, error('Intervallo di integrazione non valido.'), end
	x = linspace(a, b, n+1);
	f = feval(fun, x);
	If = ((b-a)/n)*(f(1)/2 + sum(f(2:n)) + f(n+1)/2);
end
