function x = solve_LDLt(A, x, n)
% x = solve_LDLt(A, x, n)
% La funzione risolve il sistema associato alla matrice A (di dimensione n),
% fattorizzata LDLt
% dall'algoritmo 3.6, ed al termine noto x.
% A = matrice LDLt restituita dall'algoritmo 3.6
% x = vettore dei termini noti
% n = dimensione della matrice

    x = solve_trinf_diaguni(A, x, n);
    x = solve_diag(A, x, n);
    x = solve_trisup_diaguni(A', x, n);
end
