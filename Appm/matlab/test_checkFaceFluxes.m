% Script to check data at faces on wall at z = zmin
clear
clc


%% Read mesh data
dualMesh = readMesh('dual');
fc = dualMesh.faceCenter;
idx = fc(:,3) < -2.49;
vc = dualMesh.vc;


filelist = dir('faceFluxes*.dat');
figure(2)
clf

for i = 1 : length(filelist)
filename = filelist(i).name;
fluxes = importdata(filename)';

hold on
plot(fluxes(idx,4))
hold off
title(i)
drawnow
% pause
end
grid on
return

%% Plot
figure(1)
clf
hold on
plot3(fc(idx,1), fc(idx,2), fc(idx,3), '.')
patch('Faces', dualMesh.f2v(idx,:), 'Vertices', vc, 'FaceColor', 'r', 'FaceAlpha', 0.1);
hold off
axis equal

