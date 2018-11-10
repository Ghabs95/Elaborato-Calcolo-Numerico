%Calcolo con n fissata della funzione e della sua approssimazione
n = 40;
xi = linspace(-6, 6, n+1);
fi = 1./(1+xi.^2);

x = linspace(-6,6,100000);
y1 = lagrange(xi, fi, x); %y approssimati
y2 = 1./(1+x.^2); %funzione


plot(xi,fi,'ro',x,y1,'b--',x,y2,'g');


%Calcolo della tabella con n-norma
% rst = [0 0];
% i=1;
% 
% for n=2:2:40
%     xi = linspace( -6, 6,n+1);
%     fi = 1./(1+xi.^2);
% 
%     x = linspace(-6,6,10000);
%     y1 = lagrange(xi, fi, x);
%     y2 = 1./(1+x.^2);
%     
%     rst(i,1) = n;
%     %rst(i,2) = norm(t1-y2,inf);
%     rst(i,2) = (2^(n+1))/(exp(1)*n*log(n));
%     
%     i=i+1;
% end
% 
% colNames = {'n','lebesgue'};
% tableResult = array2table(rst,...,
%     'VariableNames',colNames);
% disp(tableResult);


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