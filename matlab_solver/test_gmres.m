n = 21;
b = afun(ones(n,1));
tol = 1e-12;  maxit = 15;
p1 = 1;
p2 = 1;
x1 = gmres(@afun, b, [], [], [], [], [], [], p1, p2);

function y = afun(x, p1, p2)
n = 21;
y = [0; x(1:n-1)] + ...
[((n-1)/2:-1:0)'; (1:(n-1)/2)'].*x + ...
[x(2:n); 0];
end

