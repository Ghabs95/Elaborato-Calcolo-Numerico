function [x,i,nr] = jacobi(A, b, tol, x0, maxit)
% [x,i] = jacobi(A, b, tol, [xo, maxit])
% Restituisce la soluzione del sistema lineare Ax=b approssimata con il 
% metodo di Jacobi e il numero di iterazioni eseguite.
% Input:
%    - A: matrice utilizzata per il calcolo
%    - b: vettore dei termini noti
%    - tol: tolleranza dell' approssimazione
%    - [x0]: vettore di partenza
%    - [maxit]: numero massimo di iterazioni
% Output:
%    - x: soluzione approssimata del sistema

    % - A non deve avere elementi nulli sulla diagonale
    D = diag(A);
    if ~all(D) % all(D) ritorna vero se tutti elementi !=0
        error('La diagonale di A non deve avere elementi nulli');
    end
    n = length(b);
    if nargin <= 3
        x = rand(n,1);
    else
        x = x0;
    end
    if nargin <= 4
        maxit = 100*n*round(-log(tol));
    end
    for i = 1:maxit
        r = A * x - b;
        nr = norm(r, inf);
        if nr <= tol
            break;
        end
        r = r./D;
        x = x - r;   
    end
    if nr > tol
        warning('Raggiunto maxit.');
    end
    return
end
