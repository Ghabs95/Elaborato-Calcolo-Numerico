format longg


n=1000;
tol = 10^(-5);
[x_jacobi,i_jacobi,nr_jacobi] = jacobi(matriceSparsa(n), ones(n,1), tol, zeros(n,1));
[x_gauss_seidel,i_gauss_seidel,nr_gauss_seidel] = gauss_seidel(matriceSparsa(n), ones(n,1), tol, zeros(n,1));
semilogy(1:i_jacobi,nr_jacobi,1:i_gauss_seidel,nr_gauss_seidel);


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


function A = matriceSparsa(n)
% A = matriceSparsa(n)
% Genera la matrice quadrata sparsa nxn, con n maggiore di 10.

    %n deve essere maggiore di 10
    if n <= 10
        error('n deve essere maggiore di 10.');
    end
    
    d = ones(n,1)*4;
    A = spdiags(d,0,n,n);
    
    d = ones(n,1)*(-1);
    A = spdiags(d,1,A);
    A = spdiags(d,-1,A);
    
    if n > 10
        A = spdiags(d,10,A);
        A = spdiags(d,-10,A);
    end
end
