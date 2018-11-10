f = @(x)(x^4);
table = zeros(10,3);
x = 1;
j = 1;
while (j<=10)
    table(j, 1) = j;
    table(j, 2) = 10^(-j);
    table(j, 3) = (f(x + table(j, 2)) - f(x - table(j, 2))) / (2*table(j, 2));
    j = j + 1;
end