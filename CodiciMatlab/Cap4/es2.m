function y = newton(xi,fi,x)
% y = newton(xi,fi,x)
% Implementa il calcolo del polinomio interpolante
% di grado n in forma di Newton
% Input:
%  xi: vettore delle ascisse di interpolazione
%  fi: vettore dei valori della funzione su x
%  x: vettore dei punti in cui valutare il polinomio
%
% Output:
%  y: vettore dei valori del polinomio valutato sui punti x.
    if isempty(x), error('Vettore x vuoto!'), end
    if isempty(xi), error('Vettore xi vuoto!'), end
    if length(fi)~=length(xi), error('Vettori fi ed xi non compatibili!'), end
    dd = diffdivN(xi,fi);
    y = hornerGen(xi,dd,x);
return
end

function fi = diffdivN(xi,fi)
% f = diffdivN(xi,fi)
% Calcolo delle differenze divise per il polinomio di Newton
% Questo metodo prende in input:
%  xi: vettore delle ascisse
%  fi: vettore dei valori della funzione
%
% E restituisce:
%  fi: vettore contenente le differenze divise f[x0],f[x1,x2]...f[x0...xn]
    % Assumiamo che le ascisse siano tutte distinte.
    n = length(xi)-1;
    for j = 1:n
        for i = n+1:-1:j+1
            fi(i) = (fi(i)-fi(i-1))/(xi(i)-xi(i-j));
        end
    end
return
end

function p = hornerGen(xi,fi,x)
% p = hornerGen(xi,fi,x)
% Algoritmo di horner generalizzato
% Questo metodo prende in input:
%  fi: array dei coefficienti
%  xi: array degli zeri del polinomio
%  x: array di punti di valutazione del polinomio
% E restituisce la valutazione del polinomio
% pr-1(x)+fwr(x) nei punti x0
    n = length(xi)-1;
    p = fi(n+1)*ones(size(x));
    for k = 1:length(x) 
        for i = n:-1:1
            p(k) = p(k)*(x(k)-xi(i))+fi(i);
        end
    end
return
end
