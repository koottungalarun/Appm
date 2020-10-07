clear
clc
% close all

set(0, 'DefaultAxesFontSize', 12)

% Define file names to be read
filenames = {
    '011-explicit-coarse/011-explicit-coarse.csv', 
    '012-explicit-coarse/012-explicit-coarse.csv', 
    '013-explicit-coarse/013-explicit-coarse.csv'
    };
legnd = {'analytic', 'fine', 'medium', 'coarse'};

% Load data
for k = 1 : length(filenames)
    dataInput = importdata(filenames{k});
    idx = find(strcmp(dataInput.colheaders,'"Points:2"'));
    data(k).z = dataInput.data(:,idx);

    idx = find(strcmp(dataInput.colheaders,'"Argon density"'));
    data(k).rho = dataInput.data(:,idx);

    idx = find(strcmp(dataInput.colheaders,'"Argon pressure"'));
    data(k).p = dataInput.data(:,idx);

    idx = find(strcmp(dataInput.colheaders,'"Argon velocity:2"'));
    data(k).u = dataInput.data(:,idx);
end
dataAnalytic = load('shockTubeAnalytic.mat');

grayColor = 0.5 * [1 1 1];

figure(1)
clf
plot(dataAnalytic.x, dataAnalytic.rho, 'Color', grayColor)
hold on
% for k = 1 : length(data)
plot(data(1).z, data(1).rho, '-', 'Color', 'k')
plot(data(2).z, data(2).rho, '--', 'Color', 'k')
plot(data(3).z, data(3).rho, '-.', 'Color', 'k')
% end
hold off
grid on
ylabel('n')
ylim([0 1])
xlim([-2.5 2.5])


figure(2)
clf
plot(dataAnalytic.x, dataAnalytic.p, 'Color', grayColor)
hold on
for k = 1 : length(data)
    plot(data(k).z, data(k).p)
end
hold off
grid on
ylabel('p')
ylim([0 1])
xlim([-2.5 2.5])


figure(3)
clf
plot(dataAnalytic.x, dataAnalytic.u, 'Color', grayColor)
hold on
% for k = 1 : length(data)
plot(data(1).z, data(1).u, '-', 'Color', 'k')
plot(data(2).z, data(2).u, '-.', 'Color', 'k')
plot(data(3).z, data(3).u, '--', 'Color', 'k')
% end
hold off
grid on
ylabel('u')
xlabel('z')
ylim([0 1])
xlim([-2.5 2.5])
legend(legnd, ...
    'Location', 'NW')

h = gca;
h.YLabel.Units = 'normalized';
h.YLabel.Rotation = 0;
h.YLabel.VerticalAlignment = 'bottom';
h.YLabel.HorizontalAlignment = 'right';
h.YLabel.Position = [0 1.03];

h.XLabel.Units = 'normalized';
h.XLabel.HorizontalAlignment = 'right';
h.XLabel.VerticalAlignment = 'top';
h.XLabel.Position = [1, -0.07];

return

figure(1)
clf
plot(data(1).z, [data(1).rho data(1).p data(1).u])
hold on
set(gca, 'ColorOrderIndex', 1)
plot(data(2).z, [data(2).rho data(2).p data(2).u], '--')
hold off
grid on

dataAnalytic = load('shockTubeAnalytic.mat');
figure(2)
clf
x = dataAnalytic.x;
y(:,1) = dataAnalytic.rho;
y(:,2) = dataAnalytic.p;
y(:,3) = dataAnalytic.u;
plot(x,y)
grid on

colorGray = 0.5 * [1 1 1];
linewidth = 1;
lineWidthGray = 0.5;
% h = figure(3);
% pos = h.Position;
% oldHeight = pos(4);
% newHeight = 780;
% newBottomPos = 100;
% pos(2) = newBottomPos;
% pos(4) = newHeight;
% h.Position = pos;

figure
clf
% subplot(3,1,1)
plot(data(1).z, data(1).rho, 'LineWidth', linewidth)
hold on
plot(data(2).z, data(2).rho, 'LineWidth', linewidth)
plot(dataAnalytic.x, dataAnalytic.rho, 'Color', colorGray, 'LineWidth', linewidth)
hold off
grid on
ylabel('n')
ylim([0 1])

figure
% subplot(3,1,2)
plot(data(1).z, data(1).p, 'LineWidth', linewidth)
hold on
plot(data(2).z, data(2).p, 'LineWidth', linewidth)
plot(dataAnalytic.x, dataAnalytic.p, 'Color', colorGray, 'LineWidth', linewidth)
hold off
grid on
ylabel('p')
ylim([0 1])

figure
% subplot(3,1,3)
plot(data(1).z, data(1).u, 'LineWidth', linewidth)
hold on
plot(data(2).z, data(2).u, 'LineWidth', linewidth)
plot(dataAnalytic.x, dataAnalytic.u, 'Color', colorGray, 'LineWidth', linewidth)
hold off
ylabel('u')
grid on
ylim([0 1])
