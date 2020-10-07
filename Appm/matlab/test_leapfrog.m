% Test: Leapfrog integration

clear
clc

% Read mesh data
primal = readMesh('primal');
dual = readMesh('dual');

% Initialize state variables
nPrimalFaces = size(primal.f2e,1);
nPrimalEdges = size(primal.f2e,2);
nDualFaces = size(dual.f2e,1);
nDualEdges = size(dual.f2e,2);

bvec = zeros(nPrimalFaces,1);
evec = zeros(nPrimalEdges,1);
dvec = zeros(nDualFaces,1);
hvec = zeros(nDualEdges,1);
jvec = zeros(nDualFaces,1);

% Initialize electric current
for i = 1 : size(jvec)
    fc = dual.faceCenter(i,:);
    if vecnorm(fc(1:2)) < 0.35
        jvec(i) = 0 * dot(dual.faceNormal(i,:), [0,0,1]) * dual.faceArea(i);
    end
end

% Initialize voltage at primal vertices
nPrimalVertices = size(primal.e2v,2);
phiVec = zeros(nPrimalVertices,1);
isTerminal = logical(size(phiVec));
for i = 1 : nPrimalVertices
    pos = primal.vc(i,:);
    zPos = pos(3);
    
    isTerminal(i) = (zPos == 0 || zPos == 1) && vecnorm(pos(1:2)) < 0.35;
    phiVec(i) = zPos;
end

idx = ~isTerminal;
plot3(primal.vc(idx,1), primal.vc(idx,2), primal.vc(idx,3), 'b.')
hold on
idx = isTerminal;
plot3(primal.vc(idx,1), primal.vc(idx,2), primal.vc(idx,3), 'r.')
hold off
grid on
axis equal

G = primal.e2v;
evec = -G * phiVec;

vec = evec ./ primal.edgeLength .* primal.edgeDir;

quiver3(primal.edgeCenter(:,1), primal.edgeCenter(:,2), primal.edgeCenter(:,3), ...
    vec(:,1), vec(:,2), vec(:,3), 0)



return





%%
Cprim = primal.f2e;
Cdual = dual.f2e;

temp = dual.edgeLength(1:nPrimalFaces) ./ primal.faceArea;
hodgeOp_nu = spdiags(temp, 0, nPrimalFaces, nPrimalFaces);

temp = primal.edgeLength ./ dual.faceArea(1:nPrimalEdges);
hodgeOp_epsInv = spdiags(temp, 0, nPrimalEdges, nPrimalEdges);
clear temp


%% Time loop
dT = 0.05;
time = 0;
maxIter = 1000;
maxTime = 20;
iter = 1;
while iter <= maxIter && time <= maxTime 
    bvec = bvec - dT * Cprim * evec;
    hvec(1:size(hodgeOp_nu,1)) = hodgeOp_nu * bvec;
    dvec = dvec + dT * Cdual * hvec - dT * jvec;
    evec = hodgeOp_epsInv * dvec(1:size(hodgeOp_epsInv,2));
    time = time + dT;
    
    E = evec ./ primal.edgeLength;
    fprintf('max(abs(E)): %f \n', max(abs(E)));
    
    scale = 1;
    
    % plot primal edge values (E-field)
    figure(1)
%     subplot(2,2,1)
    pos = primal.edgePos;
    vec = evec .* primal.edgeDir ./ vecnorm(primal.edgeDir);
    quiver3(pos(:,1), pos(:,2), pos(:,3), vec(:,1), vec(:,2), vec(:,3), scale)
    axis equal
    title(sprintf('E-field, time = %f', time))
    xlabel('x')
    ylabel('y')
    zlabel('z')
%     
%     % plot dual edge values (H-field)
%     figure(2)
% %     subplot(2,2,2)
%     pos = dual.edgePos;
%     vec = hvec .* dual.edgeDir ./ vecnorm(dual.edgeDir);
%     quiver3(pos(:,1), pos(:,2), pos(:,3), vec(:,1), vec(:,2), vec(:,3))
%     axis equal
%     title(sprintf('H-field, time = %f', time))
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     
%     
%     % plot dual face values (D-field, J-field)
%     figure(3)
% %     subplot(2,2,3)
%     pos = dual.faceCenter;
%     vec = dvec .* dual.faceNormal ./ dual.faceArea;
%     quiver3(pos(:,1), pos(:,2), pos(:,3), vec(:,1), vec(:,2), vec(:,3), 'DisplayName', 'D')
%     hold on
%     vec = jvec .* dual.faceNormal ./ dual.faceArea;
%     quiver3(pos(:,1), pos(:,2), pos(:,3), vec(:,1), vec(:,2), vec(:,3), 'DisplayName', 'J')
%     hold off
%     axis equal
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
%     title(sprintf('D-/J-field, time = %f', time))
%     % legend show
%     
%     % plot primal face values (B-field)
%     figure(4)
% %     subplot(2,2,4)
%     pos = primal.faceCenter;
%     vec = bvec .* primal.faceNormal ./ primal.faceArea;
%     quiver3(pos(:,1), pos(:,2), pos(:,3), vec(:,1), vec(:,2), vec(:,3))
%     axis equal
%     title(sprintf('B-field, time = %f', time))
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
    
    drawnow
    
    iter = iter + 1;
    time = time + dT;
end






