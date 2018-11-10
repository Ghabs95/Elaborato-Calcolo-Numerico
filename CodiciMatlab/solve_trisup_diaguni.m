function x = solve_trisup_diaguni(A, x, n)
% x = trisup_diaguni(A, x, n)
% Il metodo risolve il sistema lineare associato alla matrice A di
% dimensione n triangolare superiore a diagonale unitaria e vettore dei
% termini noti x che in seguito sara' il vettore soluzione.
% A = matrice triangolare superiore a diagonale unitaria
% x = vettori dei termini noti e della soluzione
% n = dimensione della matrice


for i = n : -1 : 1
    for j = i + 1 : n
        x(i) = x(i) - A(i, j) * x(j);
    end
end