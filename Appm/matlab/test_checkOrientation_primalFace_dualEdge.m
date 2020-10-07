% Check orientation of primal face and dual edge

clear
clc

primal = readMesh('primal');
dual = readMesh('dual');

nPf = size(primal.f2e,1); % number of primal faces

idx = 1 : nPf;
pFaceN = primal.faceNormal;
dEdgeDir = dual.edgeDir;

oriented = zeros(nPf,1);
for i = idx
    A = pFaceN(i,:);
    A = A / vecnorm(A);
    B = dEdgeDir(i,:);
    B = B / vecnorm(B);
    oriented(i) = dot(A,B);
end
assert(all(oriented > 0.999))

figure(1)
clf
plotEdges(dual, 0, 'k')
hold on
idx = 1;
patch('Faces', primal.f2v(idx,:), 'Vertices', primal.vc, 'FaceColor', 'r', 'FaceAlpha', 0.1)
quiver3(primal.faceCenter(idx,1), primal.faceCenter(idx,2), primal.faceCenter(idx,3), ...
    primal.faceNormal(idx,1), primal.faceNormal(idx,2), primal.faceNormal(idx,3), 'Color', 'm', 'LineWidth', 1)

quiver3(dual.edgePos(idx,1), dual.edgePos(idx,2), dual.edgePos(idx,3), ...
    dual.edgeDir(idx,1), dual.edgeDir(idx,2), dual.edgeDir(idx,3), 'Color', 'r', 'LineWidth', 2)

hold off
grid on
axis equal
