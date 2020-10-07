% Test for DEC of APPM model, using 1st order time integration (leap frog)

clear
clc

primal = readMesh('primal');
dual = readMesh('dual');

%%
nPv = size(primal.e2v,2);
nPvt = sum(primal.vertexType == 2); % terminal vertices
nPvb = sum(primal.vertexType == 1 | primal.vertexType == 2);
nPvi = sum(primal.vertexType == 3);
assert(sum([nPvi, nPvb]) == nPv)

nPe = size(primal.e2v,1);
nPei = sum(primal.edgeType <= 1); % inner edges
nPeb = sum(primal.edgeType == 2); % boundary edges
assert(sum([nPei, nPeb]) == nPe);

nPf = size(primal.f2e,1);
nPfi = sum(primal.isFaceBoundary == 0);

%% Operators G (grad), C (curl)
G = primal.e2v;
C = primal.f2e;
Cdual = dual.f2e(1:nPe, 1:nPf);

%% Operators Meps, Mnu
values = dual.faceArea(1:nPe) ./ primal.edgeLength;
Meps = spdiags(values, 0, nPe, nPe);
Meps_inv = spdiags(1./values, 0, nPe, nPe);
values = dual.edgeLength(1:nPf) ./ primal.faceArea;
Mnu = spdiags(values, 0, nPf, nPf);

%% Operator X, GX, Q
X = vertcat(sparse(nPvi, nPvb), speye(nPvb));
GX = G*X;
Q = horzcat( vertcat(speye(nPei), sparse(nPeb,nPei)), -GX );

%% Maxwell equations with 1st order time integration (Leapfrog)
% 

dt = 1;
iter = 1;
time = dt;
iterMax = 10;
timeMax = 10;

phi1 = 0;
phi2 = @(t) 0.5 * (1 + tanh( (t - 3)));

%%
H_h = zeros(nPf,1);
E_h = zeros(nPe,1);
phiB_h = zeros(nPvb,1);

while (iter <= iterMax && time <= timeMax)
    phi_t = repmat([phi1, phi2(time)], nPvt/2,1);
    phiB_h(end+1-(1:nPvt)) = phi_t(:);
    
    E_h = E_h - -GX*phiB_h - dt * Meps_inv * Cdual * H_h;
    H_h = H_h - dt * Mnu * C * E_h;
    
    vec = E_h ./ primal.edgeLength .* primal.edgeDir;
    pos = primal.edgeCenter - 0.0 * vec;
    quiver3(pos(:,1), pos(:,2), pos(:,3), ...
        vec(:,1), vec(:,2), vec(:,3), 0)
    title(sprintf('iter = %d, time = %f', iter, time))
    drawnow
    iter = iter + 1;
    time = time + dt;
end
