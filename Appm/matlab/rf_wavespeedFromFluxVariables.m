clear
clc

syms p n u e gamma htot
syms f1 f2 f3

htot = 1/2 * u^2 + p/n* gamma/(gamma-1);
c2 = gamma * p / n;

expr1 = subs(c2, p, f2 - n*u^2);
expr2 = subs(expr1, n*u, f1);
expr3 = subs(expr2, n*u*htot, f3)
