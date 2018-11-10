xi = [0,1];
fi = [1,2];
fi1 = [2,1];
x = [0,1,2,-2,5,6,7];

hermite(xi, fi, fi1, x)

function y = hermite(xi,fi,fi1,x)
% y = hermite(xi,fi,f1i,x)
% Calcola il polinomio interpolante di Hermite
% Questo metodo prende in input:
%  xi: vettore delle ascisse di interpolazione
%  fi: vettore dei valori della funzione su x
%  fi1: vettore della derivata su x
%  x: vettore dei punti in cui valutare il polinomio
%
% E restituisce:
%  y: vettore dei valori del polinomio valutato sui punti x.
    if isempty(x), error('Vettore x vuoto!'), end
    if isempty(xi), error('Vettore xi vuoto!'), end
    if length(xi)~=length(fi) || length(xi)~=length(fi1)
        error('Vettori xi, fi e fi1 di lunghezza diversa!')
    end
    n = length(xi)-1;
    xi = reshape([xi; xi], [], 1)';
    fi = reshape([fi; fi1], [], 1)';
    dd = diffdivH(xi,fi);
    y = hornerGen(xi,dd,x);
return
end

function fi = diffdivH(xi,fi)
% fi = diffdivH(xi,fi)
% Calcola differenze divise per polinomio di Hermite
% Questo metodo prende in input:
%  xi: vettore delle ascisse
%  fi: vettore dei valori della funzione
%
% E restituisce:
%  fi: vettore contenente le differenze divise f[x0],f[x0,x0],...,f[x0,...,xn]
    n = length(xi)-1; %2*m+1
    for i = n:-2:3
        fi(i) = (fi(i)-fi(i-2))/(xi(i)-xi(i-1));
    end
    for j = 2:n
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