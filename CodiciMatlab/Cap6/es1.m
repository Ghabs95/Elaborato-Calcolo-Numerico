function A = matriceSparsa(n)
% A = matriceSparsa(n)
% Genera la matrice quadrata sparsa nxn, con n maggiore di 10.
%   
% Input:
%    - n: numero di righe/colonne della matrice quadrata sparsa
% Output:
%    - A: matrice quadrata nxn sparsa
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
