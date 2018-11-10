function [x,i,nr] = gauss_seidel(A, b, tol, x0, maxit)
% [x,i] = jacobi(A, b, tol, [xo, maxit])
% Restituisce la soluzione del sistema lineare Ax=b approssimata con il 
% metodo di Gauss-Seidel e il numero di iterazioni eseguite.
% Input:
%    - A: matrice utilizzata per il calcolo
%    - b: vettore dei termini noti
%    - tol: tolleranza dell' approssimazione
%    - [x0]: vettore di partenza
%    - [maxit]: numero massimo di iterazioni
% Output:
%    - x: soluzione approssimata del sistema
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
        err = norm(r,inf);
        nr(i) = err;
        if err<=tol
            break;
        end
        r = Msolve(A,r);
        x = x-r;
    end
    if err>tol
        warning('Non raggiunta tolleranza richiesta.');
    end
    return
end

function u = Msolve(M, r)
    % u = Msolve(M, r)
    % Restituisce la soluzione del sistema lineare Mx=r.
    u = r;
    n = length(u);
    for i = 1:n
       u(i) = u(i)/M(i,i);
       u(i+1:n) = u(i+1:n) - M(i+1:n,i)*u(i);
    end
    return
end
