% Test Rusanov flux

clear
clc

gamma = 1.4;

% Fluid state: rho, rho*u, rho*etot
primitiveLeft.rho = 1.0;
primitiveLeft.p = 1.0;
primitiveLeft.u = 0;

primitiveRight.rho = 0.125;
primitiveRight.p = 0.1;
primitiveRight.u = 0;

prim2state = @(p) [p.rho, p.rho*p.u, p.rho * (0.5 * p.u^2 + 1/(gamma-1) * p.p/p.rho)]';
qL = prim2state(primitiveLeft);
qR = prim2state(primitiveRight);

maxWavespeed = @(p) sqrt(gamma * p.p / p.rho);

sL = maxWavespeed(primitiveLeft);
sR = maxWavespeed(primitiveRight);
s = max(sL, sR);

fluxfcn = @(q) [
    q(2);
    (gamma-1) * q(3) + (3-gamma)/2 * q(2)^2/q(1);
    gamma * q(3) * q(2)/q(1) + (1-gamma)/2 * q(2)^3/q(1)^2;
    ];

fL = fluxfcn(qL);
fR = fluxfcn(qR); 
flux = 0.5 * (fL + fR) - 0.5 * s * (qR - qL);
disp(flux)

dx = 1;
dt = dx / s;

