% Test Leapfrog with Electric Potential Boundary Conditions
% 2.12.2019

clear
clc

%% Load primal and dual mesh
primalMesh = readMesh('primal');
dualMesh = readMesh('dual');

edgeType = primalMesh.edgeType;
isBoundaryEdge = edgeType == 2;
isBoundaryVertex = logical(primalMesh.isBoundaryVertex);

nPrimalEdges = size(primalMesh.edgeType,1);
nE_inner = sum(primalMesh.edgeType == 0);
nE_bn    = sum(primalMesh.edgeType == 1);
nE_notBoundary = sum(edgeType ~= 2);
nV_notBoundary = sum(primalMesh.isBoundaryVertex == 0);
nV_boundary    = sum(primalMesh.isBoundaryVertex == 1);
assert(sum(~isBoundaryEdge) == nE_notBoundary);

ndof = nE_notBoundary + nV_boundary;
nBV_offset = nE_notBoundary - nV_notBoundary;

%% Setup of operator Q
triplets = [];
for edgeIdx = 1 : nPrimalEdges
    [~,vIdxA] = find(primalMesh.e2v(edgeIdx,:) == -1);
    [~,vIdxB] = find(primalMesh.e2v(edgeIdx,:) == +1);
    if edgeType(edgeIdx) == 0
        assert(primalMesh.isBoundaryVertex(vIdxA) == 0)
        assert(primalMesh.isBoundaryVertex(vIdxB) == 0)
        triplets(end+1,:) = [edgeIdx, edgeIdx, 1];
        
    elseif edgeType(edgeIdx) == 1
        assert(sum(primalMesh.isBoundaryVertex([vIdxA, vIdxB])) == 1)
        triplets(end+1,:) = [edgeIdx, edgeIdx, 1];
        if primalMesh.isBoundaryVertex(vIdxA)
            triplets(end+1,:) = [edgeIdx, vIdxA + nBV_offset, +1];
        elseif primalMesh.isBoundaryVertex(vIdxB)
            triplets(end+1,:) = [edgeIdx, vIdxB + nBV_offset, -1];
        else
            assert(false)
        end
        
    elseif edgeType(edgeIdx) == 2
        assert(primalMesh.isBoundaryVertex(vIdxA) == 1)
        assert(primalMesh.isBoundaryVertex(vIdxB) == 1)
        triplets(end+1,:) = [edgeIdx, vIdxA + nBV_offset, +1];
        triplets(end+1,:) = [edgeIdx, vIdxB + nBV_offset, -1];
    else
        assert(false);
    end
end

assert(all(triplets(:,1) >= 1))
assert(all(triplets(:,1) <= nPrimalEdges))
assert(all(triplets(:,2) >= 1))
assert(max(triplets(:,2) == ndof))
assert(all(triplets(:,2) <= ndof))
assert(all(abs(triplets(:,3)) <= 1))

% Matrix that maps from degrees-of-freedom to primal edges: e = Q*x
Q = sparse(triplets(:,1), triplets(:,2), triplets(:,3), nPrimalEdges, ndof);
spy(Q)

%%
x = zeros(ndof,1);
pos = primalMesh.vc(isBoundaryVertex,:);

% index in x: non-boundary edges, terminal1 vertices, and terminal2 vertices
idx_EnB = (1:ndof)' <= nE_notBoundary; 
idx_V_t1 = [false(nE_notBoundary,1); vecnorm(pos(:,1:2)')' <= 0.35 & pos(:,3) <= 0.05];
idx_V_t2 = [false(nE_notBoundary,1); vecnorm(pos(:,1:2)')' <= 0.35 & pos(:,3) >= 0.95];
idx_V_nt = ~idx_EnB & ~idx_V_t1 & ~idx_V_t2; % index in x: non-terminal vertices
idx_V = idx_V_t1 | idx_V_t2 | idx_V_nt; % index in x: all boundary vertices
assert(all(size(idx_V_t1) == size(x)))
assert(all(size(idx_V_t2) == size(x)))
assert(all(sum([idx_EnB idx_V_t1 idx_V_t2 idx_V_nt],2) == 1))


phi = primalMesh.vc(:,3) + 1;
plot(primalMesh.vc(:,3), phi, '.')

G = primalMesh.e2v;
evec = -G * phi;

pos = primalMesh.edgePos;
vec = evec ./ primalMesh.edgeLength .* primalMesh.edgeDir;
quiver3(pos(:,1), pos(:,2), pos(:,3), ...
    vec(:,1), vec(:,2), vec(:,3), 0)
