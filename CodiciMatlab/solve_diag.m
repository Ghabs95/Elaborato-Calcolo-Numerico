function x = solve_diag(A, x, n)
% x = solve_diag(A, x, n)
% La funzione risolve il sistema associato alla matrice A (di dimensione n), 
% diagonale, ed al termine noto x.
% A = matrice diagonale
% x = vettore dei termini noti
% n = dimensione della matrice

   for i = 1 : n
        x(i) = x(i) / A(i, i);
   end
end