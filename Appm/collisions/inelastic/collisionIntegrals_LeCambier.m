% Calculate collision integrals from cross sections.
%
% Reference: Le, Cambier (2016).
% Modeling of inelastic collisions in a multifluid
% plasma: Ionization and recombination
% https://doi.org/10.1063/1.4953050

% Raw data of inelastic cross sections for Argon:
% Yanguas-Gil (2005): doi:10.1088/0022-3727/38/10/014
% and references therein.
%
% The ionization cross section is taken from:
% Rapp D and Englander-Golden P 1965 J. Chem. Phys. 43 1464


%% Integral for friction coefficient
% We have 
%
% $$ I_{R_0^{ion}} = \int_{x^*}^{\infty} dx_0
%        x_0^2 \exp(-x_0)
%        \zeta^{(1)} (\sqrt(\lambda x_0)) \sigma_{ion}(e), 
%        \quad e = x_0 k_B T_e 
% $$
% 
% with
% 
% $$ \zeta^{(1)} (\xi) =
%       \frac{3}{4\xi^2} (\cosh(2\xi) - \frac{1}{2\xi} \sinh(2\xi) ) 
% $$
% 
% and
% $$ \lim_{\xi \to 0} \zeta^{(1)} (\xi) = 1 $$
%
% The cross-sections are tabulated by the electron kinetic energy e
% (= impact energy) of the collision. 
%
% Note that the cross sections are 
% tabulated in units of eV, not J. Hence, we should add the conversion
% factor and find $$ e = x_0 k_B T_e (eV/J) $$.
% 
clear
clc

% Boltzmann constant
kB = 1.380649e-23; 

% unit conversion factors
eV_to_J = 1.602176634e-19; 
J_to_eV = 1 / eV_to_J;     

% ionization energy in units of electronVolt (eV) and Joule (J)
eStar_eV = 15.76; 
eStar_J = eStar_eV * eV_to_J;


%% Form function $$ \zeta^{(1)} $$
xi = linspace(0,3);
zeta = zeta1(xi);

figure(1)
plot(xi,zeta)
set(gca, 'YScale', 'log')
grid on
xlabel('\xi')
ylabel('\zeta^{(1)}(\xi)')

%% factor $$ x_0^2 * \exp(-x_0) $$
n = 200;
x0 = linspace(0, 50, n);
y = x0.^2 .* exp(-x0);

figure(1)
plot(x0, y)
xlim([0 20])
grid on
title('lin-scaled')

figure(2)
plot(x0, y)
set(gca, 'YScale', 'log')
grid on
title('log-scaled')

%% Plot cross section

% Note that the cross section data have an energy threshold at 
% ionization energy e = 15.76 eV

% n = 200;
% e_lo = 10; % lower energy bound
% e_hi = 1000; % upper energy bound
% e = logspace(log10(e_lo),log10(e_hi), n);
[~,csData] = getCrossSection([]);
e = csData(:,1);
cs = csData(:,2);

figure(1)
plot(e, cs)
grid on
title('lin-scaled')
xlabel('e / eV')
ylabel('\sigma / m^2')

figure(2)
semilogx(e, cs)
grid on
title('log-scaled')
xlabel('e / eV')
ylabel('\sigma / m^2')


%% Plot integrand without cross section

clear x0 Te lambda
fun = @(x0, Te, lambda) x0.^2 .* exp(-x0) .* zeta1(sqrt(lambda*x0));

n = 200;
x0 = logspace(-3,2,n);
lambda = 0;
Te = 0;
f = fun(x0, Te, lambda);
semilogy(x0, f)
grid on


%% Plot integral without cross section
clear x0 y
% electron temperature (Kelvin)
TeVec = linspace(300, 30e3);
TeVec = TeVec(:);
y = zeros(size(TeVec));
xStarVec = zeros(size(y));

for k = 1 : length(TeVec)
    Te = TeVec(k);
    xStar = eStar_J /(kB * Te);
    y(k) = integral(@(x0) fun(x0, Te, lambda), xStar, Inf);
    xStarVec(k) = xStar;
end
figure(1)
plot(TeVec, y)
grid on
xlabel('x_0')
ylabel('x_0^2 exp(-x_0)')



%% Plot integral

clear x0 Te lambda y xStarVec
fun = @(x0, Te, lambda) ...
    x0.^2 .* exp(-x0) ...
    .* zeta1(sqrt(lambda*x0)) ...
    .* getCrossSection(x0 * kB * Te * J_to_eV);

n = 50;
lambdaVec = [0 logspace(-2, 1, n)];
% lambdaVec = [0 0.01 0.1 1 5];
% lambdaVec = lambdaVec(:);
TeVec = TeVec(:);


y = zeros(size(TeVec, 1), size(lambdaVec, 1));
xStarVec = zeros(size(y));

for m = 1 : length(lambdaVec)
    lambda = lambdaVec(m);
    for k = 1 : length(TeVec)
        Te = TeVec(k);
        xStar = eStar_J / (kB * Te);
        x_up = 10000;
        y(k,m) = integral(@(x0) fun(x0, Te, lambda), xStar, x_up);
        xStarVec(k,m) = xStar;
    end
    lgnd{m} = sprintf('\\lambda = %.2e', lambda);
end


if length(lambdaVec) <= 5
figure(1)
semilogy(TeVec, y)
xlabel('T_e / K')
ylabel('I / m^2')
legend(lgnd, 'Location', 'SE')
grid on
ylim([1e-40 1e-15])

figure(2)
plot(TeVec, y)
xlim(20e3 + 5e3*[-1 1])
ylim([1e-23 1e-22])
legend(lgnd, 'Location', 'SE')
grid on
xlabel('T / K')
ylabel('I / m^2')

% It is useful to consider the data in terms of 1/T.
% See, e.g., Annaloro (?)
figure(3)
plot(1./TeVec, log10(y))
xlabel('(1/T) / (K^{-1})')
ylabel('log_{10} (I / m^2)')
legend(lgnd, 'Location', 'NE')
grid on
end

% Plot surface data: $$ I = I(T_e, \lambda) $$
figure(4)
surf(lambdaVec, 1./TeVec, log10(y), 'EdgeAlpha', 0)
xlabel('\lambda')
ylabel('T^{-1} / (K^{-1})')
zlabel('I / m^2')

figure(5)
contourf(lambdaVec, 1./TeVec, log10(y));
xlabel('\lambda')
ylabel('T^{-1} / K^{-1}')
colorbar
title('log10(I / m^2)')

figure(6)
imagesc(lambdaVec, 1./TeVec, log10(y))
xlabel('\lambda')
ylabel('T^{-1} / K^{-1}')
colorbar

%% Write data to file
M = [NaN lambdaVec;
    TeVec, y];
filename = 'e_Ar_inelasticCrossSectionIntegrals R0ion.csv';
fid = fopen(filename, 'w');
assert(fid >= 3);
fprintf(fid, '# Inelastic cross sections integrals R0ion \n');
fprintf(fid, '# cf. eq. B1a in Le and Cambiert (2016), https://doi.org/10.1063/1.4953050 \n')
fprintf(fid, '# \n')
fprintf(fid, '# Data source: Rapp D and Englander-Golden P 1965 J. Chem. Phys. 43 1464 \n');
fprintf(fid, '# Format: \n');
fprintf(fid, '#   first row: lambda (ratio of electron kinetic energy to thermal energy \n')
fprintf(fid, '#   first col: electron temperature (K) \n');
fprintf(fid, '#   data: cross section integral (m^2) for R0ion \n')
fprintf(fid, '# \n');
fclose(fid);
writematrix(M, filename, 'WriteMode', 'append');

