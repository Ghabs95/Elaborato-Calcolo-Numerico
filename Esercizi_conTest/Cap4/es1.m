xi = [0,1,2,3];
fi = [1,1,2,2];
x = [1,2,3,4,5,6,7,0];

lagrange(xi, fi, x)

function y = lagrange(xi, fi, x)
% y = lagrange(xi,fi,x)
% Implementa il calcolo del polinomio interpolante
% di grado n in forma di Lagrange
% Questo metodo prende in input:
%  xi: vettore delle ascisse di interpolazione
%  fi: vettore dei valori della funzione su x
%  x: vettore dei punti in cui valutare il polinomio
%
% E restituisce:
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