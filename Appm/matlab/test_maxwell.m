clear
clc

primal = readMesh('primal');
dual = readMesh('dual');

iter = 200;
filename = sprintf('appm-%d.h5', iter);

H = h5read(filename, '/H')';
H_h = h5read(filename, '/hvec');
pos = dual.edgeCenter;
vec = H;

figure(1)
quiver3(pos(:,1), pos(:,2), pos(:,3), ...
    vec(:,1), vec(:,2), vec(:,3))

nPf = size(primal.f2e,1);
values = primal.faceArea ./ dual.edgeLength(1:nPf);
M_mu = spdiags(values, 0, nPf, nPf);
B_h = M_mu * H_h(1:nPf);

pos = primal.faceCenter;
vec = repmat(B_h ./ primal.faceArea, 1, 3) .* primal.faceNormal;

figure(2)
quiver3(pos(:,1), pos(:,2), pos(:,3), ...
    vec(:,1), vec(:,2), vec(:,3))

