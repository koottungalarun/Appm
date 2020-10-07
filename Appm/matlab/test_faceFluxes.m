% Test: check sum of face fluxes. It seems that numerical noise leads to
% unstable results.

% It doesn't help. Check first that the primal mesh is as accurate as
% possible!

% error('');

clear
% clc

%% Load mesh data
dualMesh = readMesh('dual');
cell_sumMomentumFlux = h5read('appm-1.h5', '/Neutrals-sumMomentumFlux')';
faceArea = dualMesh.faceArea;
faceFluxes = h5read('appm-1.h5', '/Neutrals-momentumFlux')';
cc = dualMesh.cellCenter';
nCells = size(dualMesh.c2f,1);
nFaces = size(dualMesh.c2f,2);
cellVolume = dualMesh.cellVolume;

%% Select a cell in fluid domain
summed = zeros(nCells, 3);
for idx = 1 : nCells
    if dualMesh.cellType(idx) ~= 1
        continue
    end
    % Get indices of adjacient faces
    [~,cellFaceIdx,~] = find(dualMesh.c2f(idx,:));
    
    c = 0;
    for i = cellFaceIdx
        ff = faceFluxes(i,:);
%         ff = single(faceFluxes(i,:));
        incidence = full(dualMesh.c2f(idx, i));
        inp = incidence * ff * faceArea(i);

%         Ordinary sum
        summed(idx,:) = summed(idx,:) + inp;

        % Kahan's summation algorithm: improves on numerical cancellation
%         y = inp - c;
%         t = summed(idx,:) + y;
%         c = (t - summed(idx,:)) - y;
%         summed(idx,:) = t;
        
    end
    summed(idx,:) = summed(idx,:) / cellVolume(idx);
%     break
end
plot(summed, '.')
all(summed(:) == 0)
all(cell_sumMomentumFlux(:) == 0)
