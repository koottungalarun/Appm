
clear
clc

syms gamma
syms q1 q2 q3
syms f1 f2 f3

assume(f1, 'real')
assume(f2, 'real')
assume(f3, 'real')
assume(gamma > 1)

eqn1 = f1 == q2;
eqn2 = f2 == (gamma-1) * q3 + (3-gamma)/2 * q2^2/q1;
eqn3 = f3 == gamma * q2 * q3 / q1 - (gamma-1)/2 * q2^3 / q1^2;


%% Substitute eqn1 (q2 -> f1) into eqn2 and eqn3
eqn2 = subs(eqn2, rhs(eqn1), lhs(eqn1));
eqn3 = subs(eqn3, rhs(eqn1), lhs(eqn1));

%% Isolate equations for q1 and q3
% Isolate q1 from eqn2 to lhs
eqn2 = isolate(eqn2, q1);

% Isolate q3 from eqn3 to lhs
eqn3 = isolate(eqn3, q3);

%% Substitute eqn3 into eqn2 and solve for q1
eqn2 = subs(eqn2, lhs(eqn3), rhs(eqn3));
pretty(eqn2)
eqn2 = isolate(eqn2, q1);
pretty(eqn2)

%% Substitute this expression for q1 into eqn3, and solve for q3
eqn3 = subs(eqn3, lhs(eqn2), rhs(eqn2));
pretty(eqn3)

%% Substitute with actual values
subs(eqn3, [f1, f2, f3, gamma], [1, 2, 4, 7/5])
subs(eqn2, [f1, f2, f3, gamma], [1, 2, 4, 7/5])

syms p csq
eqn_p = p == (gamma - 1) * (q3 - 1/2 * f1^2/q1);
eqn_csq = csq == gamma * p * q1/f1;

eqn_csq = subs(eqn_csq, lhs(eqn_p), rhs(eqn_p));
eqn_csq = subs(eqn_csq, lhs(eqn2), rhs(eqn2));
eqn_csq = subs(eqn_csq, lhs(eqn3), rhs(eqn3));

assume(gamma > 1)

pretty(eqn_csq)
% pretty(simplify(eqn_csq))

assume(f1 ~= 0)
isolate(simplify(eqn_csq), csq)
pretty(isolate(simplify(eqn_csq), csq))

% or:
% csq = gamma * (-2 + 1/(2*(gamma-1)) * f2/f3 
%             * (1 + sqrt(1 - 2*f1*f3/(f2^2) * (gamma^2 - 1) / gamma^2))

return
subs(eqn_csq, [f1 f2 f3 gamma], [1 2 4 7/5])


