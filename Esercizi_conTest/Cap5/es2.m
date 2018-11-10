f = @(x) sin(x);
If_vettoriale = simpcomp( 2, 0, pi, f)

function If = simpcomp(n, a, b, fun)
% If = simpcomp(n, a, b, fun)
% Calcola l'integrale della funzione, nell'intervallo prescelto,
% usando la formula di Simpson composita.
% Input:
 % 	n: intero positivo pari, indica num intervalli in [a,b]
 % 	a: estremo sinistro
 % 	b: estremo destro
 % 	fun: funzione integranda
% Output:
 % 	If: approssimazione dell'integrale definito della funzione
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