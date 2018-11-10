format longg

tol = 10^(-5);
for i=1:10
    n=i*100;
    rst(i,1) = n;
    [rst(i,3) rst(i,2)] = potenze(matriceSparsa(n), tol, ones(n,1));
end


colNames = {'n','numero_iterazioni','stima_autovalore'};
tableResult = array2table(rst,'VariableNames',colNames);
disp(tableResult);


% A = ...
% [2 -1 0 0 0 0 0 0 0 0;
% -1 2 -1 0 0 0 0 0 0 0;
% 0 -1 2 -1 0 0 0 0 0 0;
% 0 0 -1 2 -1 0 0 0 0 0;
% 0 0 0 -1 2 -1 0 0 0 0;
% 0 0 0 0 -1 2 -1 0 0 0;
% 0 0 0 0 0 -1 2 -1 0 0;
% 0 0 0 0 0 0 -1 2 -1 0;
% 0 0 0 0 0 0 0 -1 2 -1;
% 0 0 0 0 0 0 0 0 -1 2];

% tol = 10^(-5);

% l = potenze(A, tol)


function [lambda, i] = potenze(A, tol, x0, maxit)
% [lambda, i] = metodoPotenze(A, tol, [x0, maxit])
% Restituisce l'autovalore dominante della matrice A e il numero di 
% iterazioni necessarie per calcolarlo.
%   
% Input:
%    - A: matrice utilizzata per il calcolo
%    - tol: tolleranza dell' approssimazione
%    - [x0]: vettore di partenza
%    - [maxit]: numero massimo di iterazioni
% Output:
%    - lambda: matrice quadrata nxn sparsa
%    - i: numero di iterazioni
    [m,n] = size(A);
    if m ~= n 
        error('La matrice deve essere quadrata.'); 
    end
    if x0(:)==0
        error('Il vettore x0 non può avere esclusivamente elementi nulli.');
    end % if da eliminare per calcolo con matrice
    
    if nargin <= 2
        x = rand(n, 1); 
    else
        x = x0;
    end
    if nargin <= 3
        maxit = 100*2*round(-log(tol));
    end
    x = x0; % da eliminare per calcolo con matrice
    x = x / norm(x);
    lambda = inf;
    for i=1:maxit
        lambda0 = lambda;
        v = A * x;
        lambda = x' * v;
        err = abs(lambda - lambda0);
        if err <= tol
            break
        end
        x = v/norm(v);
    end
    if err > tol
        warning("Raggiunto maxit.");
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