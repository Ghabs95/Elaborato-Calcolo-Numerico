function y = lagrange(xi, fi, x)
% y = lagrange(xi,fi,x)
% Implementa il calcolo del polinomio interpolante
% di grado n in forma di Lagrange
% Input:
%  xi: vettore delle ascisse di interpolazione
%  fi: vettore dei valori della funzione su x
%  x: vettore dei punti in cui valutare il polinomio
%
% Output:
%  y: vettore dei valori del polinomio valutato sui punti x
   if isempty(x), error("Vettore x vuoto!"), end
   if length(xi)~=length(fi), error("Diversa lunghezza dei vettori xi e fi!"), end
    n=length(xi)-1;
    m=length(x);
    y=zeros(size(x));
    for i=1:m
        for j=1:n+1
            p=1;
            for k=1:n+1
                if j~=k
                    p=p*(x(i)-xi(k))/(xi(j)-xi(k));
                end
            end
            y(i)=y(i)+fi(j)*p;
        end
    end
return
end
