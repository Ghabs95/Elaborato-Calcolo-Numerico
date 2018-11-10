function x = solve_LU(A, x, n)
% x = solve_LU(A, x, n)
% La funzione risolve il sistema associato alla matrice A (di dimensione n),
% fattorizzata LU
% dall'algoritmo 3.7, ed al termine noto x.
% A = matrice LU restituita dall'algoritmo 3.7
% x = vettore dei termini noti
% n = dimensione della matrice

    x = solve_trinf_diaguni(A, x, n);
    x = solve_trisup(A, x, n);
end