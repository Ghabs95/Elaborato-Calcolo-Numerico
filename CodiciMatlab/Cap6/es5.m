n = 1000;
tol = 10^(-5);
[x_jacobi,i_jacobi,nr_jacobi] = jacobi(matriceSparsa(n), ones(n,1), tol, zeros(n,1));
[x_gauss_seidel,i_gauss_seidel,nr_gauss_seidel] = gauss_seidel(matriceSparsa(n), ones(n,1), tol, zeros(n,1));
semilogy(1:i_jacobi,nr_jacobi,1:i_gauss_seidel,nr_gauss_seidel);