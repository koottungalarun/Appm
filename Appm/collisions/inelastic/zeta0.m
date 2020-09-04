function zeta = zeta0(xi)

zeta = zeros(size(xi));

tol = 1e-2;
idx = xi > tol;
x = xi(idx);
zeta(idx) = sinh(2*x) ./ (2*x);

% Taylor approximation
% syms x
% f = sinh(2*x) / (2*x)
% T = taylor(f, 'Order', 6)
% T = (2*x^4)/15 + (2*x^2)/3 + 1

x = xi(~idx);
zeta(~idx) = 1 + 2/3 * x.^2 + 2/15 * x.^4;

end