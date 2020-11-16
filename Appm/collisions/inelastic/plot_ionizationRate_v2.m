clear
clc

%% Physics constants
eV_to_J = 1.60218e-19;   % conversion from eV to J
h = 6.626e-34;           % Planck constant, Js
kB = 1.38064852e-23;     % Boltzmann constant, J/K
me = 9.10938e-31;        % electron mass, kg
E_ion = 15.76 * eV_to_J; % ionization energy of argon (in eV, to J)


%% Ionization rate as computed by Le & Cambier (2016) with data form Rapp (1965)
dataIon = readmatrix('I_Gion.csv');
TeVecIon = dataIon(2:end,1);
I_Gion = dataIon(2:end,2);
ve = sqrt(8*kB*TeVecIon/(pi*me));

ki = ve .* I_Gion;

%% (Direct) Ionization rate by approximate formula (Benilov & Naidis, 1998)
ci = 18e-22 / eV_to_J;
Ry = 13.6 * eV_to_J;
E = 4.11 * eV_to_J;
dE = E_ion - E;
kdir = ci * (8/pi*kB/me * TeVecIon).^0.5 .* (E_ion + 2 * kB * TeVecIon) .* exp(-E_ion ./ (kB * TeVecIon));

%% Recombination rate as computed by Le & Cambier (2016) and data from Rapp (1965)
dataRec = readmatrix('I_Grec.csv');
TeVec = dataRec(2:end,1);
assert(all(TeVec == TeVecIon))
I_Grec = dataRec(2:end,2);

xStar = E_ion ./ (kB * TeVec);
lambda = h ./ sqrt(2*pi*me*kB*TeVec);
ve = sqrt(8*kB*TeVec/(pi*me));

kr = 1/2 * 1/6 * lambda.^3 .* ve .* exp(xStar) .* I_Grec;

%% Plot data
figure(3)
clf

Tscale = 1e3;
yyaxis left
ph(2) = semilogy(TeVecIon / Tscale, ki, 'DisplayName', 'k_{ion}');
hold on
ph(1) = plot(TeVecIon / Tscale, kdir, 'DisplayName', 'k_{dir}');
hold off
ylim([1e-35 1e-15])
ylabel('k_{ion}, k_{dir} / (m^3 s^{-1})')

ax = gca;
ax.YLabel.Units = 'normalized';
ax.YLabel.Position = [-0.05 1.03];
ax.YLabel.Rotation = 0;
ax.YLabel.HorizontalAlignment = 'left';
ax.YLabel.VerticalAlignment = 'bottom';

x = 0.28 + [0.1 0.0];
y = 0.85 * [1 1];
t_h = annotation('textarrow', x, y, 'String', 'k_{ion}, k_{dir}')

yyaxis right
ph(3) = semilogy(TeVec / Tscale, kr, 'DisplayName', 'k_{rec}', 'LineStyle', '-.');
grid on
xlabel('T / kK')
ylabel('k_{rec} / (m^6 s^{-1})')

x = 0.75 + [0.0 0.1];
y = 0.3 * [1 1];
t_h = annotation('textarrow', x, y, 'String', 'k_{rec}');


ax = gca;
legend(ph);
ax.Legend.Visible = 'on';
ax.Legend.Location = 'East';

ax = gca;
ax.YLabel.Units = 'normalized';
ax.YLabel.Position = [1.05 1.03];
ax.YLabel.Rotation = 0;
ax.YLabel.HorizontalAlignment = 'right';
ax.YLabel.VerticalAlignment = 'bottom';

ax.XLabel.VerticalAlignment = 'top';
ax.XLabel.HorizontalAlignment = 'right';
ax.XLabel.Units = 'normalized';
ax.XLabel.Position = [1 -0.05];


isPrintFigure = true;
if isPrintFigure
    print('-deps2', 'Ar-reactionRates-ionization-recombination.eps')
end
return
%% Check with Saha equation
saha = 2./lambda.^3 * 6 .* exp(-E_ion ./ (kB * TeVec));
figure(4)
clf
semilogy(TeVec, saha, 'DisplayName', 'Saha eq. (analytic)')
hold on
semilogy(TeVec, ki ./ kr, 'DisplayName', 'ki / kr (= Saha)')
hold off
legend show
ylim([1e10 1e30])
grid on
title('Saha')