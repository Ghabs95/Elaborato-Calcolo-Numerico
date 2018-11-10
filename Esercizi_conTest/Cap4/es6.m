xi = ceby(2, 1, 5)

function xi = ceby(n, a, b)
% xi = ceby(n, a, b)
% Implementa il calcolo delle ascisse di Chebyshev
% per il polinomio di grado n su [a,b]
% Questo metodo prende in input:
%  n: grado del polinomio interpolante
%  a: estremo sinistro dell'intervallo
%  b: estremo destro dell'intervallo
%
% E restituisce:
%  xi: ascisse di Chebyshev
    if n<=0
        error('Inserire numero maggiore di 0!')
    end
    xi = cos((2*[0:n]+1)*pi/(2*n+2));
    xi = ((a+b)+(b-a)*xi)/2;
    return
end
