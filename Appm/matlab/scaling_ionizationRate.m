% Validity check for scaling of ionization rate and recombination rate

clear
clc

% Physics constants
kB = 1.38e-23; % J / K, Boltzmann constant
c0 = 3e8; % m/s, speed of light
eps0 = 8.85e-12;  % F/m
q = 1.602e-19; % C, charge scale
oneBarPressure = 10135; % Pa, unit conversion bar->Pa
oneElectronVolt_Kelvin = 11604;
eV_to_J = q;

xbar = 1e-2; % m, length scale
Tbar = 1 * oneElectronVolt_Kelvin; % K, temperature scale
pbar = 1e0 * oneBarPressure; % Pa, pressure scale
nbar = pbar / (kB * Tbar); % m^-3, number density scale by ideal gas law

E_i_dim = 15.75; % eV, ionization energy
E_i = E_i_dim * eV_to_J / (kB * Tbar);



% AP-parameter
lambda = sqrt(eps0 * kB * Tbar / (q^2 * nbar * xbar^2));

% Time scale
tbar = 1/c0 * xbar / lambda;

% Reaction rate scale
kbar = 1/(nbar * tbar);

data = importdata('collisions/inelastic/e-Ar.dat');

Te_dim = data(:,1);
ki_dim = data(:,2);
kr_dim = data(:,3);

Te = Te_dim / Tbar;
ki = ki_dim / kbar;
kr = kr_dim / (kbar / nbar);

figure(1)
semilogy(Te, ki, 'DisplayName', 'k_i')
hold on
semilogy(Te, kr, 'DisplayName', 'k_r')
hold off
grid on
legend show


