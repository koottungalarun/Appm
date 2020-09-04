% Get collision integrals

clear all
clc

%% Consistent cross section data and its derivative
% Check that data and its derivative are consistent
n = 1e4;
e = [0 logspace(-3, log10(1e3), n)];

% get data
cs = getCrossSection(e);
ds = getCrossSectionDrv(e);

% find ionization energy in derivative of cross section data;
% i.e., where derivative sufficiently gets positive
idx = find(ds > 1e-23, 1);
E_ion = e(idx);
assert(E_ion >= 15.75 & E_ion <= 15.77)

figure(1)
nSubPlots = 2;
subplot(nSubPlots,1,1)
plot(e, cs, 'DisplayName', '\sigma(e)')
legend('Location', 'NE')
grid on
xlabel('e / eV')
ylabel('\sigma / m^2')

subplot(nSubPlots,1,2)
plot(e, ds, 'DisplayName', 'd\sigma(e)/de')
legend('Location', 'NE')
grid on
xlabel('e / eV')
ylabel('(d\sigma/de) / (m^2 / eV)')

figure(2)
plot(e, cumtrapz(e,ds), 'Color', 'r', 'DisplayName', '\int d\sigma/de de')
hold on
plot(e, cs, 'Color', 'k', 'LineStyle', '--', 'DisplayName', '\sigma')
hold off
legend('Location', 'NE')
xlabel('e / eV')
ylabel('\sigma / m^2')
grid on

%% Physics constants

% Boltzmann constant
kB = 1.380649e-23; 

% unit conversion factors
eV_to_J = 1.602176634e-19; 
J_to_eV = 1 / eV_to_J;     

% ionization energy in units of electronVolt (eV) and Joule (J)
eStar_eV = 15.76; 
eStar_J = eStar_eV * eV_to_J;


%% Parameters
TeVec = (300 : 100 : 40e3)';
lambdaVec = [0 0.01 0.1 0.3 1 3 10]; %logspace(-2, 1, 50)];

%% Integrands
fun_Gion = @(x0, Te) ...
    x0 .* exp(-x0) ...
    .* getCrossSection(x0 * kB * Te * J_to_eV);

fun_R0ion = @(x0, Te, lambda) ...
    x0.^2 .* exp(-x0) .* zeta1(sqrt(lambda*x0)) ...
    .* getCrossSection(x0 * kB * Te * J_to_eV);

fun_J00ion = @(x0, Te, lambda) ...
    x0.^2 .* exp(-x0) .* zeta0(sqrt(lambda*x0)) ...
    .* getCrossSection(x0 * kB * Te * J_to_eV);

fun_Grec = @(x0,v,xStar,Te,lambda) ...
    x0 .* exp(-x0) ...
    .* zeta0(sqrt(lambda * (x0 - v))) .* zeta0(sqrt(lambda * (v - xStar))) ...
    .* getCrossSectionDrv(x0 * kB * Te * J_to_eV);

fun_R1rec = @(x0,v,xStar,Te,lambda) ...
    x0 .* exp(-x0) .* (x0 - v) ...
    .* zeta1(sqrt(lambda*(x0-v))) .* zeta0(sqrt(lambda*(v-xStar))) ...
    .* getCrossSectionDrv(x0*kB*Te*J_to_eV);
    
fun_R2rec = @(x0,v,xStar,Te,lambda) ...
    x0 .* exp(-x0) .* (v - xStar) ...
    .* zeta0(sqrt(lambda*(x0-v))) .* zeta1(sqrt(lambda*(v-xStar))) ...
    .* getCrossSectionDrv(x0*kB*Te*J_to_eV);

fun_J11rec = @(x0,v,xStar,Te,lambda) ...
    x0 .* exp(-x0) .* (x0 - v) ...
    .* zeta0(sqrt(lambda*(x0-v))) .* zeta0(sqrt(lambda*(v-xStar))) ...
    .* getCrossSectionDrv(x0*kB*Te*J_to_eV);

fun_J22rec = @(x0,v,xStar,Te,lambda) ...
    x0 .* exp(-x0) .* (v - xStar) ...
    .* zeta0(sqrt(lambda*(x0-v))) .* zeta0(sqrt(lambda*(v-xStar))) ...
    .* getCrossSectionDrv(x0*kB*Te*J_to_eV);

fun_J12rec = @(x0,v,xStar,Te,lambda) ...
    x0 .* exp(-x0) .* (x0 - v) .* (v - xStar) ...
    .* zeta1(sqrt(lambda*(x0-v))) .* zeta1(sqrt(lambda*(v-xStar))) ...
    .* getCrossSectionDrv(x0*kB*Te*J_to_eV);



%% Compute integrals

x_upBound = 10000;
I_Gion = zeros(size(TeVec));
I_R0ion = zeros([length(TeVec), length(lambdaVec)]);
I_Grec = zeros(size(I_R0ion));
I_R1rec = zeros(size(I_R0ion));
I_R2rec = zeros(size(I_R0ion));
I_J00ion = zeros(size(I_R0ion));
I_J11rec = zeros(size(I_R0ion));
I_J22rec = zeros(size(I_R0ion));
I_J12rec = zeros(size(I_R0ion));

timer1 = tic;
for m = 1 : length(lambdaVec)
    lambda = lambdaVec(m);
    fprintf('m = %d / %d', m, length(lambdaVec));
    timer2 = tic;
    for k = 1 : length(TeVec)
        Te = TeVec(k);
        xStar = eStar_J / (kB * Te);
        if m == 1
            I_Gion(k) = integral(@(x0) fun_Gion(x0, Te), xStar, x_upBound);
        end
        I_R0ion(k,m) = integral(@(x0) fun_R0ion(x0, Te, lambda), xStar, x_upBound);
        I_J00ion(k,m) = integral(@(x0) fun_J00ion(x0, Te, lambda), xStar, x_upBound);
        I_Grec(k,m) = integral2(@(x0,v) fun_Grec(x0,v,xStar,Te,lambda), xStar, x_upBound, xStar, @(x0) x0);
        I_R1rec(k,m) = integral2(@(x0,v) fun_R1rec(x0,v,xStar,Te,lambda), xStar, x_upBound, xStar, @(x0) x0);
        I_R2rec(k,m) = integral2(@(x0,v) fun_R2rec(x0,v,xStar,Te,lambda), xStar, x_upBound, xStar, @(x0) x0);
        I_J11rec(k,m) = integral2(@(x0,v) fun_J11rec(x0,v,xStar,Te,lambda), xStar, x_upBound, xStar, @(x0) x0);
        I_J22rec(k,m) = integral2(@(x0,v) fun_J22rec(x0,v,xStar,Te,lambda), xStar, x_upBound, xStar, @(x0) x0);
        I_J12rec(k,m) = integral2(@(x0,v) fun_J12rec(x0,v,xStar,Te,lambda), xStar, x_upBound, xStar, @(x0) x0);
    end
    t_m = toc(timer2);
    fprintf('\t time = %f \n', t_m);
end
t_total = toc(timer1);
sprintf('Total time: %f\n', t_total);


%% Plot data
% Lidx = 1 : 10 : length(lambdaVec);
Lidx = 1 : length(lambdaVec);
clear lgnd
for i = 1 : length(Lidx)
    lgnd{i} = sprintf('\\lambda = %.2f', lambdaVec(Lidx(i)));
end

data{1} = I_R0ion;
data{2} = I_J00ion;
data{3} = I_R1rec;
data{4} = I_R2rec;
data{5} = I_J11rec;
data{6} = I_J22rec;
data{7} = I_J12rec;
data{8} = I_Gion;
data{9} = I_Grec;

titleStr = [
    "R0ion", 
    "J00ion",
    "R1rec",
    "R2rec",
    "J11rec",
    "J22rec",
    "J12rec", 
    "Gion",
    "Grec"
    ];
assert(length(data) == length(titleStr));

for i = 1 : min(7,length(data))
    figure(i)
    clf
    semilogy(1./TeVec, data{i})
    grid on
    ylabel('I / m^2')
    xlabel('(T/K)^{-1}')
    legend(lgnd)
    title(titleStr(i))
    ylim([1e-50 1e-10])
end

marker = [
    "o",
    "s",
    "d",
    "x",
    "*",
    "^",
    "v"
    ];

figure(8)
clf
idx = 1;
for i = 1 : min(7,length(data))
    hold on
    semilogy(1./TeVec, data{i}(:,idx), ...
        'Marker', marker(i), ...
        'MarkerIndices', [1 11 22 33 44 55], ...
        'DisplayName', titleStr(i))
    hold off
end
set(gca, 'YScale', 'log')
legend show
grid on
ylim([1e-50 1e-10])
title(sprintf('\\lambda = %f', lambdaVec(idx)))

%%
figure(9)
clf
semilogy(1./TeVec, I_Gion, 'Color', 'k', 'DisplayName', 'I \Gamma^{ion}')
hold on
for i = 1 : size(I_Grec,2)
    legName = sprintf('%.2e', lambdaVec(i));
%     legName = '';
    semilogy(1./TeVec, I_Grec(:,i), 'DisplayName', strcat('I \Gamma^{rec} \lambda=', legName))
end
hold off
set(gca, 'YScale', 'log')
legend('Location', 'best')
ylim([1e-50 1e-20])
grid on

%% Write data to file

for i = 1 : length(data)
    filename = strcat('I_', titleStr(i), '.csv');
    % Write header section
    fid = fopen(filename, 'w');
    fprintf(fid, '# \n'); % TODO add header that describes the data set
    fclose(fid);
    
    M = [NaN, lambdaVec(1:size(data{i},2));
        TeVec, data{i}];
    writematrix(M, filename, 'WriteMode', 'append');
end
