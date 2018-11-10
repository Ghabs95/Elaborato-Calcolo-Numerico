function y = spline3( xi, fi, x, tipo )
% y = spline3( xi, fi, x, tipo )
% Calcola e valuta i valori della spline cubica naturale o not-a-knot 
% (rispettivamente per tipo uguale a 0 o 1) delle coppie di dati assegnati.
% Input:
%    xi: ascisse di interpolazione
%    fi: valori della funzione, valutati nelle ascisse xi
%    x: punti di valutazione della spline 
%    tipo: 
%        - nat identifica la spline naturale
%        - nak identifica la spline not-a-knot
%        - un valore diverso restituisce errore
% Output:
%    y: valutazione dei punti x calcolati sulla spline richiesta
    if tipo == "nat"
        mi = sistemaSplineNaturale(xi, fi);
        y = valutaSplineNaturale(xi, fi, mi, x);
    elseif tipo == "nak"
        mi = sistemaSplineNotAKnot(xi, fi);        
        y = valutaSplineNotAKnot(xi, fi, mi, x);
    else
        error('Il tipo indicato non è corretto.');
    end
    return
end

function fin = differenzeDivise(x, f)
% Input:
% - f: funzione;
% - x: vettore delle ascisse.
% Output:
% - fin: l'ultima differenza divisa calcolata, corrispondente a quella che
%        include tutte le ascisse.
    n = length(x);
    for i = 1:n-1
        for j = n:-1:i+1
            f(j)=(f(j)-f(j-1))/(x(j)-x(j-i));
        end
    end
    fin = f(end);
end

function m = sistemaSplineNaturale( xi , fi )
% Input:
%    xi: ascisse di interpolazione
%    fi: valori della funzione, valutati nelle ascisse xi
% Output:
%    m: coefficienti necessari a calcolare l'espressione della spline
%    cubica

    % calcolo delle phi, zi e diff_div (lunghi n-1)
    % n:= grado del polinomio interpolante
    n = length(xi)-1;
    phi = zeros(n-2, 1);
    zi = zeros(n-2, 1);
    dd = zeros(n-1, 1);
    hi = xi(2) - xi(1);
    hi1 = xi(3) - xi(2);
    phi(1) = hi/(hi+hi1);
    zi(1) = hi1/(hi+hi1);
    dd(1) = differenzeDivise(xi(1:3), fi(1:3));
    for i = 2:n-2 
        hi = xi(i) - xi(i-1);
        hi1 = xi(i+1) - xi(i);
        phi(i) = hi/(hi+hi1);
        zi(i) = hi1/(hi+hi1);
        dd(i) = differenzeDivise(xi(i-1:i+1), fi(i-1:i+1));
    end
    dd(i+1) = differenzeDivise(xi(i:i+2), fi(i:i+2));
    % calcolo delle li ed ui che completano il sistema lineare della spline
    % naturale
    u(1) = 2;
    for i = 2:n-1 
        l(i) = phi(i-1)/u(i-1);
        u(i) = 2-l(i)*zi(i-1);
    end
    % risoluzione del sistema lineare LU per trovare gli m
    % 1) L*y = 6*b
    % il vettore dd delle diffDiv viene sovrascritto
    dd(1) = 6*dd(1);
    dd(2:n-1) = 6*dd(2:n-1)-l(2:n-1)*dd(1:n-2);
    % 2) U*m = y
    m(n-1) = dd(n-1)/u(n-1);
    m(n-2:-1:1) = (dd(n-2:-1:1)-zi(n-2:-1:1) * m(n-1:-1:2))/u(n-2:-1:1);
end

function y = valutaSplineNaturale(xi, fi, mi, x)
    % Input:
    %    xi: ascisse di interpolazione
    %    fi: valori della funzione, valutati nelle ascisse xi
    %    mi: coefficienti necessari a calcolare l'espressione della spline
    %    x: punti in cui valutare la spline
    % Output:
    %    y: valutazioni nei punti x della spline.
    
    if xi(1)>x(1) || xi(end)<x(end)
        error('I punti di valutazione non rientrano nel dominio della spline.');
    end
    if length(xi)~=length(fi)
        error('I vettori xi e fi devono avere la stessa lunghezza.');
    end
    % raccolgo tutti i sottoinsiemi di punti da valutare con le relative
    % functions, i punti x vegono ordinati
    mi = [0  mi  0];
    sort(x);
    y = zeros(length(x),1);    
    lastIndex = 1;
    k = 2;
    for j = 1 : length(x)
        if x(j) <= xi(k-1)
            j = j + 1;
        else
            if j ~= lastIndex
                % calcolo la spline relativa al k-1° intervallo
                hi = xi(k) - xi(k-1);
                ri = fi(k-1) - ((hi^2)/6) * mi(k-1);
                qi = (fi(k) - fi(k-1))/hi - (hi/6) * ...
                    (mi(k) - mi(k-1));
                spline = @(x) (((x - xi(k-1)).^3) .* mi(k) + ((xi(k) - x).^3) .* mi(k-1))./(6*hi)...
                  + qi .* (x - xi(k-1)) + ri;
                % calcolo le valutazioni della spline
                y(lastIndex : j-1) = spline(x(lastIndex : j-1));
                lastIndex = j;
            end
            k = k+1;
        end
    end
    % valutazione degli ultimi punti
    if j ~= lastIndex
        hi = xi(end) - xi(end-1);
        ri = fi(end-1) - ((hi^2)/6) * mi(end-1);
        qi = (fi(end) - fi(end-1))/hi - (hi/6) * (-mi(end-1));
        spline = @(x) (((x - xi(end-1)).^3) .* mi(end) + ((xi(end) - x).^3) .* mi(end-1))./(6*hi)...
                  + qi .* (x - xi(end-1)) + ri;
        y(lastIndex:j-1) = spline(x(lastIndex:j-1));
    end
    
end

function m = sistemaSplineNotAKnot(xi, fi)
% Calcola i coefficienti m da applicare all'espressione della spline cubica
% not-a-knot.
% Input:
%    xi: ascisse di interpolazione
%    fi: valori della funzione, valutati nelle ascisse xi
% Output:
%    m: coefficienti necessari a calcolare l'espressione della spline
%    not-a-knot
    
    % calcolo delle phi, zi e diff_div
    % n:= grado del polinomio interpolante
    n = length(xi)-1;
    hi = xi(2) - xi(1);
    hi1 = xi(3) - xi(2);
    phi(1) = hi/(hi+hi1);
    zi(1) = hi1/(hi+hi1);
    diffDiv(1) = differenzeDivise(xi(1:3), fi(1:3));
    for i = 2:n-1 
        hi = xi(i) - xi(i-1);
        hi1 = xi(i+1) - xi(i);
        phi(i) = hi/(hi+hi1);
        zi(i) = hi1/(hi+hi1);
        diffDiv(i) = differenzeDivise(xi(i-1:i+1), fi(i-1:i+1));
    end
    diffDiv(i+1) = differenzeDivise(xi(i:i+2), fi(i:i+2));
    % espando diffDiv a vettore di lunghezza n+1
    diffDiv = [6*diffDiv(1), 6*diffDiv, 6*diffDiv(end)];
    % i = 1
    u(1) = 1;
    % i = 2
    w(1) = 0;
    l(1) = phi(1) / u(1);
    u(2) = 2 - phi(1);
    % i = 3
    w(2) = zi(1) - phi(1);
    l(2) = phi(2) / u(2);
    u(3) = 2 - l(2)*w(2);
    % i = 4:n-1
    for i = 4:n-1
        w(i-1) = zi(i-2);
        l(i-1) = phi(i-1) / u(i-1);
        u(i) = 2 - l(i-1)*w(i-1);
    end
    % i = n
    w(n-1) = zi(n-2);
    l(n) = (phi(n-1) - zi(n-1)) / u(n-1);
    u(n) = 2 - zi(n-1) - l(n)*w(n-1);
    % i = n+1
    w(n) = zi(n-1);
    l(n+1) = 0;
    u(n+1) = 1;
    % 1) Ly = 6*diffDiv
    y(1) = diffDiv(1);
    for i = 2:n+1
        y(i) = diffDiv(i) - l(i-1)*y(i-1);
    end
    % 2) Um = y
    m = zeros(n+1, 1); 
    m(n+1) = y(n+1) / u(n+1);
    for i = n:-1:1
        m(i) = (y(i) - w(i) * m(i+1)) / u(i);
    end
    m(1) = m(1) - m(2) - m(3);
    m(end)= m(n+1)-m(n)-m(n-1);
    return
end

function y = valutaSplineNotAKnot(xi, fi, mi, x)
% Valuta i punti della spline cubica not-a-knot nei punti assegnati x
% Input:
%    xi: ascisse di interpolazione
%    fi: valori della funzione, valutati nelle ascisse xi
%    mi: coefficienti necessari a calcolare l'espressione della spline
%    x: punti in cui valutare la spline
% Output:
%    y: valutazioni nei punti x della spline.
    
    if xi(1)>x(1) || xi(end)<x(end)
        error('I punti di valutazione non rientrano nel dominio della spline.');
    end
    if length(xi)~=length(fi)
        error('I vettori xi e fi devono avere la stessa lunghezza.');
    end
    % raccolgo tutti i sottoinsiemi di punti da valutare con le relative
    % functions, i punti x vegono ordinati
    sort(x);
    y = zeros(length(x),1);    
    lastIndex = 1;
    k = 2;
    for j = 1 : length(x)
        if x(j) <= xi(k-1)
            j = j + 1;
        else
            if j ~= lastIndex
                % calcolo la spline relativa al k-1° intervallo
                hi = xi(k) - xi(k-1);
                ri = fi(k-1) - ((hi^2)/6) * mi(k-1);
                qi = (fi(k) - fi(k-1))/hi - (hi/6) * ...
                    (mi(k) - mi(k-1));
                spline = @(x) (((x - xi(k-1)).^3) .* mi(k) + ((xi(k) - x).^3) .* mi(k-1))./(6*hi)...
                  + qi .* (x - xi(k-1)) + ri;
                % calcolo le valutazioni della spline
                y(lastIndex : j-1) = spline(x(lastIndex : j-1));
                lastIndex = j;
            end
            k = k+1;
        end
    end
    % valutazione degli ultimi punti
    if j ~= lastIndex
        hi = xi(end) - xi(end-1);
        ri = fi(end-1) - ((hi^2)/6) * mi(end-1);
        qi = (fi(end) - fi(end-1))/hi - (hi/6) * (-mi(end-1));
        spline = @(x) (((x - xi(end-1)).^3) .* mi(end) + ((xi(end) - x).^3) .* mi(end-1))./(6*hi)...
                  + qi .* (x - xi(end-1)) + ri;
        y(lastIndex:j-1) = spline(x(lastIndex:j-1));
    end
end
