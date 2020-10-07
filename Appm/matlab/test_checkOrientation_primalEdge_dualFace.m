% Compare orientation of primal edge and dual face

clear
clc

primal = readMesh('primal');
dual = readMesh('dual');


pEdgeDir = primal.edgeDir(1,:);
pEdgeC = primal.edgeCenter(1,:);

dFaceN = dual.faceNormal(1,:);
dFaceC = dual.faceCenter(1,:);

% number of primal edges
pNe = size(primal.e2v,1);

for i = 1 : pNe
    A = primal.edgeDir(i,:);
    A = A / vecnorm(A);
    B = dual.faceNormal(i,:);
    B = B / vecnorm(B);
    oriented(i) = dot(A, B);
end

figure(1)
clf
plotEdges(primal, 0, 'k')
idx = 2;
title(sprintf('idx = %d', idx))
hold on
quiver3(primal.edgePos(idx,1), primal.edgePos(idx,2), primal.edgePos(idx,3), ...
    primal.edgeDir(idx,1), primal.edgeDir(idx,2), primal.edgeDir(idx,3))
hold off

patch('Faces', dual.f2v(idx,:), 'Vertices', dual.vc, 'FaceColor', 'g', 'FaceAlpha', 0.1)
hold on
quiver3(dual.faceCenter(idx,1), dual.faceCenter(idx,2), dual.faceCenter(idx,3), ...
    dual.faceNormal(idx,1), dual.faceNormal(idx,2), dual.faceNormal(idx,3))
hold off

set(gca, 'View', [14 20])
grid on
axis equal

figure(2)
plot(oriented, '.')