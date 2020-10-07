clear
clc

% Read data of explicit case
data = importdata('Testcase-Euler-Shocktube-Explicit.csv');
idx(1) = find(strcmp(data.colheaders, """Points:2"""));
idx(2) = find(strcmp(data.colheaders, """Neutrals number density"""));
idx(3) = find(strcmp(data.colheaders, """Neutrals pressure"""));
idx(4) = find(strcmp(data.colheaders, """Neutrals velocity:2"""));

explicitData.z = data.data(:,idx(1));
explicitData.n = data.data(:,idx(2));
explicitData.p = data.data(:,idx(3));
explicitData.u = data.data(:,idx(4));
explicitData.data = data.data(:,idx);

% Read data of implicit case
data = importdata('Testcase-Euler-Shocktube-Implicit.csv');
idx(1) = find(strcmp(data.colheaders, """Points:2"""));
idx(2) = find(strcmp(data.colheaders, """Neutrals number density"""));
idx(3) = find(strcmp(data.colheaders, """Neutrals pressure"""));
idx(4) = find(strcmp(data.colheaders, """Neutrals velocity:2"""));

implicitData.z = data.data(:,idx(1));
implicitData.n = data.data(:,idx(2));
implicitData.p = data.data(:,idx(3));
implicitData.u = data.data(:,idx(4));
implicitData.data = data.data(:,idx);

figure(1)
clf
% set(gca, 'ColorOrderIndex', 1)
p(1:3) = plot(implicitData.z, implicitData.data(:,2:end), 'LineWidth', 2, 'LineStyle', '-');
hold on
set(gca, 'ColorOrderIndex', 1)
p(4:6) = plot(explicitData.z, explicitData.data(:,2:end), '--');
hold off
grid on

ylim([0 1.2])

xlabel('z')
legnd = ["Number density", "Pressure", "Velocity"];
legend(legnd)
