function x = solve_trinf_diaguni(A, x, n)
% x = solve_trinf_diaguni(A, x, n)
% Il metodo risolve il sistema lineare associato alla matrice A di
% dimensione n triangolare inferiore a diagonale unitaria e vettore dei
% termini noti x che in seguito sara' il vettore soluzione.
% A = matrice triangolare inferiore a diagonale unitaria
% x = vettori dei termini noti e della soluzione
% n = dimensione della matrice

for i = 1 : n
    for j = 1 : i-1
        x(i) = x(i) - A(i, j) * x(j);
    end
end