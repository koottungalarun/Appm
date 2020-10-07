% Plot data of cell at position x0 

clear
clc

%% find cell with index idx0 that is closest to position x0
x0 = [0 0 0]; % find cell with center closest to this position
dualMesh = readMesh('dual');
x = dualMesh.cellCenter - x0;
dist = vecnorm(x')';
[d,idx0] = min(dist); % idx0 is the cell index

%% get list of output files and sort by their timestamps
filenameList = dir('appm-*.h5');
[~,idx] = sort([filenameList.datenum]);
filenameList = filenameList(idx);

% check if time values in *.h5 files are sorted (increasing)
clear time
time = zeros(size(filenameList));
for i = 1 : length(filenameList)
    time(i) = h5read(filenameList(i).name, '/time');
end
assert(all(diff(time) > 0 ))
clear time

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
iterMax = length(filenameList);
for idx = 1 : iterMax
%     i = i + 100;
%     filename = sprintf('appm-%05d.h5', i);
    filename = filenameList(idx).name;
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
plot(time, 1 - n(:,3), '-.')
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


figure(2)
for i = 1 : 3
    subplot(3,1,i)
    plot(time, u(:,i))
    grid on
    legend(fluidNames(i))
    if i == 1
        title('Velocity')
    end
end

figure(3)
for i = 1 : 3
    subplot(3,1,i)
    plot(time, n(:,i))
    grid on
    legend(fluidNames(i))
    if i == 1
        title('Number density')
    end
end


