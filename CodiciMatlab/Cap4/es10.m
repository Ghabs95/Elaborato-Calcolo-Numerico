format longg


misurazioni = ...
[1 2.9;
1 3.1;
2 6.9;
2 7.1;
3 12.9;
3 13.1;
4 20.9;
4 21.1;
5 30.9;
5 31.1];


n=2;


A = vander(misurazioni(:,1));

A = A(:,8:10)

A = fliplr(A);
misurazioni(:,2);
xx = sistemaQR(fattorizzazioneQR(A),misurazioni(:,2))

poly = polyfit(misurazioni(:,1),misurazioni(:,2),2)


x = 1:6; 
f = polyval([1,1,2],x);

plot(misurazioni(:,1),misurazioni(:,2),'o',x,f,'-')



%codice della fatt QR

% Input:
% -A: matrice prodotta dalla function del precedente esercizio;
% -b: vettore dei termini noti.
% Output:
% -b: soluzioni del sistema lineare sovradeterminato.
function [b] = sistemaQR(A,b)
    [m,n] = size(A);
    Qtrasp=eye(m); %ID
    for i=1:n
        Qtrasp = [eye(i-1) zeros(i-1, m-i+1); zeros(i-1, m-i+1)' (eye(m-i+1)-(2/norm([1; A(i+1:m, i)], 2)^2)*([1; A(i+1:m, i)]*[1 A(i+1:m, i)']))]*Qtrasp;
    end
    b = triangSup(triu(A(1:n, :)), Qtrasp(1:n,:)*b);
end


function [b] = triangSup(A, b)
    for i=length(A):-1:1
        for j=i+1:length(A)
            b(i)=b(i)-A(i,j)*b(j);
        end
        if A(i,i)==0
            error('Matrice singolare')
        else
            b(i)=b(i)/A(i,i);
        end
    end
end

% Input:
% -A: Matrice 
% Output:
% -b: matrice A riscritta con la parte significativa di R e la parte 
% significativa dei vettori di Householder normalizzati con prima componente unitaria
function [A] = fattorizzazioneQR(A)
    [m,n]=size(A);
    for i=1:n
    alpha = norm(A(i:m, i), 2);
    if alpha==0
        error('La matrice non ha rank massimo');
    end
    if(A(i,i))>=0
        alpha = -alpha;
    end
    v1 = A(i,i)-alpha;
    A(i,i) = alpha;
    A(i+1:m, i) = A(i+1:m, i)/v1;
    beta = -v1/alpha;
    A(i:m, i+1:n) = A(i:m, i+1:n) -(beta*[1; A(i+1:m, i)])*([1 A(i+1:m, i)']*A(i:m, i+1:n));
    end
end