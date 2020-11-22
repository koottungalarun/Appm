% Evaluate magnetic field

% foldername = 'data/validation/emag_case2_medium_lambdaSq_1e-1';
% mkdir(foldername)
% save(strcat(foldername, '/matlab.mat'))
% copyfile('*.txt', foldername)


%%
clear
clc

dualMesh = readMesh('dual');

x0 = [0 0 0];
x = dualMesh.cellCenter - x0;
dist = vecnorm(x')';
[~,idx0] = min(dist); % idx0 is the cell index

x1 = [1 0 0];
x = dualMesh.cellCenter - x1;
dist = vecnorm(x')';
[~,idx1] = min(dist); % idx1 is the cell index

%%
filenameList = dir('appm-*.h5');
% [~,idx] = sort([filenameList.datenum]);
% filenameList = filenameList(idx);

clear time
time = zeros(size(filenameList));
for i = 1 : length(filenameList)
    time(i) = h5read(filenameList(i).name, '/time');
end
[time,idx] = sort(time);
filenameList = filenameList(idx);

if ~any(diff(time) <= 0)
    time2 = zeros(size(filenameList));
    for i = 1 : length(filenameList)
        time2(i) = h5read(filenameList(i).name, '/time');
    end
else 
    time2 = time;
end
% check if time values in *.h5 files are sorted (increasing)
assert(all(diff(time2) > 0 ))

filename = 'appm-00000.h5';
assert(isfile(filename));
info = h5info(filename);

%%
i = 0;
iterMax = length(filenameList);
% iterMax = min(100, iterMax);
% warning('iterMax is limited')
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
    for j = 1 %: length(fluidNames)
       
        if j == 1
            datasetName = '/Ecc';
            temp = h5read(filename, datasetName);
            E(1:3,idx) = temp(:,idx0);
            temp = h5read(filename, '/Bcc');
            B(1:3,idx) = temp(:,idx1);
            I_tot(1:3, idx) = h5read(filename, '/speciesTotalCurrent');
        end
    end
end



%% find cell with index idx1 that is closest to position x0
figure(1)
plot(time, E(3,:))
grid on
xlabel('t')
ylabel('E_z')


figure(2)
plot(time, B(2,:))
grid on
title('B_y at (1,0,0)')
xlabel('t')
ylabel('B_y')
