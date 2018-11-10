%Calcolo con n fissata della funzione e della sua approssimazione
n = 30;
xi = ceby(n, -6, 6);

fi = 1./(1+xi.^2); %funzione di Runge

x = linspace(-6,6,1000000);
y1 = lagrange(xi, fi, x); %y approssimati
y2 = 1./(1+x.^2); %funzione

plot(xi,fi,'ro',x,y1,'b--',x,y2,'g');

Calcolo della tabella con n-norma
rst = [0 0];
i=1;

for n=2:2:40
    xi = ceby(n, -6, 6);
    fi = 1./(1+xi.^2);

    x = linspace(-6,6,10000);
    y1 = lagrange(xi, fi, x);
    y2 = 1./(1+x.^2);
    
    
    rst(i,1) = n;
    rst(i,2) = norm(y1-y2,inf);
    
    %Crescita della costante di Lebesgue con le ascisse di Chebyshev 
    rst(i,3) = (2/pi)*log(n);
    
    i=i+1;
end

colNames = {'n','norm','Lebesgue'};
tableResult = array2table(rst,'VariableNames',colNames);
disp(tableResult);


function xi = ceby(n, a, b)
% xi = ceby(n, a, b)
% Implementa il calcolo delle ascisse di Chebyshev
% per il polinomio di grado n su [a,b]
    if n<=0
        error('Inserire numero maggiore di 0!')
    end
    xi = cos((2*[0:n]+1)*pi/(2*n+2));
    xi = ((a+b)+(b-a)*xi)/2;
    return
end

function y = lagrange(xi, fi, x)
% y = lagrange(xi,fi,x)
% Implementa il calcolo del polinomio interpolante
% di grado n in forma di Lagrange
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