% Plot mesh from indicence map

clear
clc

%% Load primal mesh data
% primal.vc = importdata('primal-vertices.dat');
% primal.f2e = readSparseMatrix('primal-f2e.dat');
% primal.e2v = readSparseMatrix('primal-e2v.dat');
% primal.f2v = importdata('primal-f2v.dat');
% primal.f2v(primal.f2v(:) < 0) = nan;
% primal.f2v = primal.f2v + 1;
% primal.c2f = readSparseMatrix('primal-c2f.dat');
primal = readMesh('primal');

assert(all(abs(primal.f2e(:)) <= 1));
assert(all(abs(primal.e2v(:)) <= 1));
assert(size(primal.f2e,2) <= size(primal.e2v, 1));
assert(size(primal.e2v,2) <= size(primal.vc,  1));
assert(size(primal.c2f,2) <= size(primal.f2e, 1));

figure(1)
clf
patch('Faces', primal.f2v, 'Vertices', primal.vc, 'FaceColor', 'r', 'FaceAlpha', 0.0)
axis equal
grid on

if any(vecnorm(primal.vc') > 1.05)
    phi = linspace(0,2*pi);
    z = exp(1i*phi);
    hold on
    plot(real(z), imag(z), 'LineWidth', 2, 'Color', 'k')
    hold off
end

showEdgeLabels = false;
color = 'r';
plotEdges(primal, showEdgeLabels, color);


%% Load dual mesh data
% clf
showDualVertexLabels = true;
showDualEdgeLabels = false;
showDualFaceLabels = true;

dual = readMesh('dual');
if isempty(dual.vc)
    fprintf('Dual mesh is empty')
    return
end
if size(dual.vc,1) > 200
    showDualVertexLabels = false;
    showDualEdgeLabels = false;
    showDualFaceLabels = false;
end

% plot vertices

hold on
plot3(dual.vc(:,1), dual.vc(:,2), dual.vc(:,3), 'ko')
if showDualVertexLabels
    text(dual.vc(:,1), dual.vc(:,2), dual.vc(:,3), num2cell( 1:size(dual.vc,1) ), 'Color', 'k')
end
hold off

% plotEdges(dual, showDualEdgeLabels, 'g')

% plot faces
patch('Faces', dual.f2v, 'Vertices', dual.vc, 'FaceColor', 'g', 'FaceAlpha', 0.1)
if showDualFaceLabels
    text(dual.faceCenter(:,1), dual.faceCenter(:,2), dual.faceCenter(:,3), num2cell( 1:size(dual.faceCenter,1) ), 'Color', 'g')
end

axis equal
grid on
shg


