function zeta = zeta1(xi)
% zeta = 3/((2x)^2) * (cosh(2*x) - 1/(2*x) * sinh(2*x))

zeta = zeros(size(xi));

tol = 1e-2;

idx = xi > tol;
x1 = xi(~idx);
zeta(~idx) = (2*x1.^4)/35 + (2*x1.^2)/5 + 1;

x2 = xi(idx);
zeta(idx) = 3./4. * x2.^(-2) .* (cosh(2*x2) - 0.5 * x2.^(-1) .* sinh(2*x2));

%% Taylor expansion
% syms x
% f = 3/(2*x)^2 * (cosh(2*x) - 1/(2*x)*sinh(2*x));
%
% T = taylor(f, 'Order', 4)
% T = (2*x^2)/5 + 1
%
% T = taylor(f, 'Order', 6)
% T = (2*x^4)/35 + (2*x^2)/5 + 1


end