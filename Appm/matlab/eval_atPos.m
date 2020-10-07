% Plot data of cell at position x0 

clear
clc

%% find cell with index idx0 that is closest to position x0
x0 = [0 0 0]; % find cell with center closest to this position
dualMesh = readMesh('dual');
x = dualMesh.cellCenter - x0;
dist = vecnorm(x')';
[d,idx0] = min(dist); % idx0 is the cell index

%% get list of output files
filenameList = dir('appm-*.h5');
iterMax = length(filenameList);

%%
fluidNames = [
    "Electron"
    "Ion"
    "Argon"
];

datasetNames = [
    "stateE"
    "pressure"
    "temperature"
    "numberDensity"
];

filename = 'appm-00000.h5';
assert(isfile(filename));
info = h5info(filename);

i = 0;
for idx = 1 : iterMax
    i = i + 100;
    filename = sprintf('appm-%05d.h5', i);
%     filename = filenameList(idx).name;
    if ~isfile(filename)
        warning("File not found: %s", filename)
        break
    end
    assert(isfile(filename));
    time(idx) = h5read(filename, '/time');
    for j = 1 : length(fluidNames)
        % Read scalar data
        start = [idx0 1];
        count = [1 Inf];

        datasetName = sprintf('/%s-%s', fluidNames(j), 'stateE');
        stateE(idx,j) = h5read(filename, datasetName, start, count);
        
        datasetName = sprintf('/%s-%s', fluidNames(j), 'pressure');
        p(idx,j) = h5read(filename, datasetName, start, count);
        
        datasetName = sprintf('/%s-%s', fluidNames(j), 'numberDensity');
        n(idx,j) = h5read(filename, datasetName, start, count);
        
        datasetName = sprintf('/%s-%s', fluidNames(j), 'temperature');
        T(idx,j) = h5read(filename, datasetName, start, count);

        % Vector data
        vecData.start = [3 idx0];
        vecData.count = [1 1];
        datasetName = sprintf('/%s-%s', fluidNames(j), 'velocity');
        u(idx,j) = h5read(filename, datasetName, vecData.start, vecData.count);
    end
end

%%
close all
pos = [50 50 1400 800]; % left bottom width height
figure('Position', pos);
rows = 1;
cols = 3;
idx = 0;

idx = idx + 1;
subplot(rows, cols, idx)
plot(time, n(:,1), '-')
hold on
plot(time, n(:,2), '--')
plot(time, n(:,3), '-.')
hold off
grid on
xlabel('t')
ylabel('n')
legend(fluidNames)
title('Number density')

idx = idx + 1;
subplot(rows,cols,idx)
plot(time, u)
grid on
xlabel('t')
ylabel('u_z')
title('Velocity')
legend(fluidNames)

idx = idx + 1;
subplot(rows,cols,idx)
plot(time, T)
grid on
xlabel('t')
ylabel('T')
title('Temperature')
legend(fluidNames)



% semilogy(time, stateE)
% grid on
% xlabel('t')
% ylabel('stateE')
% legend(fluidNames)
% 
% % subplot(3,1,2)
% figure(2)
% plot(time, p)
% grid on
% xlabel('t')
% ylabel('p')
% legend(fluidNames)


