clear 
clc

primalMesh = readMesh('primal');
dualMesh = readMesh('dual');

%%
lambda2 = 1e-2;
dt = 0.1;
iterMax = 100;
timeMax = 7;

nPv = size(primalMesh.e2v,2);
nPvb = sum(primalMesh.isBoundaryVertex);
nPvt = sum(primalMesh.vertexType == 2);

nPe = size(primalMesh.e2v,1);
nPei = sum(primalMesh.edgeType < 2);

nPf = size(primalMesh.f2e,1);
nPfi = sum(~primalMesh.isFaceBoundary);

nDf = size(dualMesh.f2e,1);

idx_edgeZdir = false(nPe,1);
for i = 1 : nPe
    edgeDir = primalMesh.edgeDir(i,:);
    idx_edgeZdir(i) = vecnorm( cross(edgeDir, [0 0 1]) ) < 1e-4;
end

B_h = zeros(nPf,1);

J_h = zeros(nDf,1);
idx_CurrentZdir = false(size(J_h));
for i = 1 : nPe %length(J_h)
    fn = dualMesh.faceNormal(i,:);
    idx_CurrentZdir(i) = vecnorm( cross(fn,[0 0 1])) < 1e-4;
end

t0 = 1;
tscale = 0.2;
currentFcn1 = @(t) 0.5 * (1 + tanh( (t - t0) / tscale ) );
currentFcn2 = @(t) 0.5 * (1 + tanh(-(t - 3)  / tscale ) );
% currentFcn = @(t) currentFcn1(t) .* currentFcn2(t);
currentFcn = @(t) exp(-((t-2)/0.5).^2) + -exp(-((t-5)/0.5).^2);


figure(1)
time = linspace(0,10);
plot(time, currentFcn(time))
grid on
clear time

t0 = 1;
tscale = 0.2;
voltageBCfunc = @(t) 0.5 * (1 + tanh( (t - t0) / tscale) );

%% Discrete operators (div, grad, curl)

C = primalMesh.f2e;
C00 = C(1:nPfi, 1:nPei);

% Q = sparse(nPe, nPei + nPvb);
X = vertcat( sparse(nPv - nPvb, nPvb), speye(nPvb, nPvb));
G = primalMesh.e2v;
Q = horzcat( speye(nPe, nPei), -G*X);
assert(size(Q,1) == nPe);
assert(size(Q,2) == nPei + nPvb);

%% Material laws
temp = dualMesh.faceArea(1:nPe) ./ primalMesh.edgeLength;
Meps = spdiags(temp, 0, nPe, nPe);

temp = dualMesh.edgeLength(1:nPfi) ./ dualMesh.faceArea(1:nPfi);
Mnu00 = spdiags(temp, 0, nPfi, nPfi);

%% System of equations
M1 = lambda2 * Q' * Meps * Q;

% P = sparse(nPfi, nPe);
% P(1:nPfi, 1:nPei) = C00;
P = horzcat(C00, sparse(nPfi, nPvb));
M2 = P' * Mnu00 * P;

A = M1 + dt^2 * M2;
rhs = 0;

%% 
x = zeros(size(Q,2), 1); % x(m+1)
x_m = x; % x(m)
x_mm1 = x; % x(m-1)

% Index vectors for Dirichlet conditions and Free indices
idx_d = false(size(x,1),1);
idx_f = false(size(x,1),1);

isCurrentSourceActive = 1;

if isCurrentSourceActive
    idx_f(1 : end - nPvt/2) = true;
else
    idx_f(1 : end - nPvt) = true;
end

idx_d = ~idx_f;

A_f = A(idx_f, idx_f);
A_d = A(:,idx_d);

iter = 0;
time = 0;
timeVec = [];

Erms_zdir = [];
while iter < iterMax && time < timeMax
    
    % Update for next timestep
    x_mm1 = x_m;
    x_m = x;
    dt_mm1 = dt;
    dt_ratio = dt / dt_mm1;
    time = time + dt;
    iter = iter + 1;
    timeVec(end+1) = time;
    
    J_h_mm1 = J_h;
    if isCurrentSourceActive
        J_h(idx_CurrentZdir) = currentFcn(time);
    end
    deltaJ = J_h - J_h_mm1;

    phiA = voltageBCfunc(time) * ones(nPvt/2,1);
    phiB = 0 * ones(nPvt/2,1);
    
    if isCurrentSourceActive
        x_d = phiB;
    else
        x_d = [phiA;
            phiB];    
    end
    
    rhs = zeros(size(x));
    rhs = rhs + (1 + dt_ratio) * M1 * x_m - dt_ratio * M1 * x_mm1;
    
    if isCurrentSourceActive
        temp = dt * Q' * deltaJ(1:nPe);
%         temp((nPei + 1) : end) = 0;
        rhs = rhs - temp;
    end

    rhs = rhs - A_d * x_d;
    rhs_f = rhs(idx_f);
    x_f = A_f \ rhs_f;
    x = [x_f; x_d];
    
    E_h = Q*x;
    
    B_h = B_h - dt * C * E_h;
    Brms(iter) = sqrt(sum(B_h.^2) / length(B_h));
    
    temp = E_h ./ primalMesh.edgeLength;
    temp = temp(idx_edgeZdir);
    Erms_zdir(iter) = sqrt(sum(temp.^2) / length(temp));

    figure(1)
    nsubfigs = 4;
    clf
    subplot(nsubfigs,1,1)
    plot(x,'.')
    grid on
    ylabel('x')
    title(sprintf('iter = %d', iter))
    
    subplot(nsubfigs,1,2)
    plot(E_h ./ primalMesh.edgeLength, '.')
    grid on
    ylabel('E')
    
    subplot(nsubfigs,1,3)
    plot(timeVec, Erms_zdir)
    ylabel('E rms')
    grid on
    
    subplot(nsubfigs,1,4)
    plot(timeVec, currentFcn(timeVec))
    ylabel('J(t)')
    grid on
    
    figure(4)
    subplot(3,1,1)
    plot(B_h ./ primalMesh.faceArea, '.')
    grid on
    ylabel('B')
    
    subplot(3,1,2)
    plot(timeVec, Brms)
    grid on
    ylabel('B rms')
    
    subplot(3,1,3)
    plot(deltaJ, '.')
    ylabel('dJ')
    grid on
    
    
    pos = primalMesh.edgePos;
    temp = E_h ./ primalMesh.edgeLength;
    vec = primalMesh.edgeDir .* repmat(temp, 1, 3);

    figure(2)
    clf
    quiver3(pos(:,1), pos(:,2), pos(:,3), ...
        vec(:,1), vec(:,2), vec(:,3))
    axis equal
    grid on
    
    figure(3)
    temp = B_h ./ primalMesh.faceArea;
    vec = primalMesh.faceNormal .* repmat(temp, 1, 3);
    pos = primalMesh.faceCenter;
    quiver3(pos(:,1), pos(:,2), pos(:,3), ...
        vec(:,1), vec(:,2), vec(:,3), 0)
    title('B field')
    grid on
    axis equal
    
    drawnow
    
    pause
    
end


%%
