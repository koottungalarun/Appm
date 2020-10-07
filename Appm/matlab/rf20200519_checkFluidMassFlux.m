clear
clc

faceIdx = 23488 + 1;
dualMesh = readMesh('dual');


%%
filename = 'appm-500.h5';
info = h5info(filename);
data = h5read(filename, '/Electrons-massFlux');

faceIdx = 23488 + 1;



%%

for iter = 1 : 500
    filename = sprintf('appm-%d.h5', iter);
    data = h5read(filename, '/Electrons-massFlux');
    massflux(iter) = data(faceIdx);
end

plot(diff(massflux))
