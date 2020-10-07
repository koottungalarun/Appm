% Test for Implicit Euler method to solve Maxwell's equations

clear
clc
primal = readMesh('primal');
dual = readMesh('dual');

%%
lambdaSquare = 1e-0;
lambda = sqrt(lambdaSquare);

% discrete curl operator
C = primal.f2e;

% discrete grad operator
G = primal.e2v;

nPvb = sum(primal.vertexType <= 2);
nPvt = sum(primal.vertexType == 2);
nPvi = sum(primal.vertexType == 3);
nPei = sum(primal.edgeType <= 1);
nPeb = sum(primal.edgeType == 2);
nPe = size(G,1);
nPfi = sum(primal.isFaceBoundary == 0);

% discrete inclusion operator for boundary vertices
X = [sparse(nPvi, nPvb);
    speye(nPvb)];

GX = G * X;
Q = [vertcat(speye(nPei), sparse(nPeb, nPei)), -GX];

Meps = spdiags(dual.faceArea(1:nPe) ./ primal.edgeLength(1:nPe), 0, nPe, nPe);
Mnu  = spdiags(primal.faceArea(1:nPfi) ./ dual.edgeLength(1:nPfi), 0, nPfi, nPfi);


zUnit = [0 0 1];
a = primal.edgeDir(1:nPei,:) ./ primal.edgeLength(1:nPei);
b = repmat(zUnit,nPei,1);
temp = dot(a',b')';

idx = abs(temp) > 0.01;
Ei = zeros(nPei,1);
Ei(idx) = primal.edgeLength(idx);

idx = primal.vertexType <= 2;
phiB = 1 - primal.vc(idx,3);
% state vector of inner edges and boundary vertices
x_h = [Ei; phiB];

E_h = Q * x_h;

C00 = C(1:nPfi, 1:nPe);

% temp = C00 * Q;
% temp = Mnu * temp;
A2 = Q' * C00' * Mnu * C00 * Q;
A1 = lambda^2 * Q' * Meps * Q;

x_m = zeros(size(x_h));
x_mm1 = x_m;
x = x_m;

nFree = size(Ei);
nD = size(phiB);

iterMax = 5000;
timeMax = 50;
iter = 0;
time = 0;
rms = [];
while iter < iterMax && time < timeMax
    % update state vectors
    x_mm1 = x_m;
    x_m = x;
    
    % define Dirichlet boundary conditions
    x_D = zeros(size(phiB));
    voltageB = 1;
    x_D(end/2 : end) = voltageB;
    
    % Define timestep
    dt = 1e-0;
    % Define system of equations
    A = A1 + dt^2 * A2;
    
    A_free = A(1:nFree, 1:nFree);
    A_D = A(:,end-nD+1 : end);
    
    rhs = lambdaSquare * Q' * Meps * Q * (2*x_m - x_mm1);
    rhs = rhs - A_D * x_D;
    
    rhs_free = rhs(1:nFree);
    x_free = A_free \ rhs_free;
    
    x = [x_free; x_D];
    
    iter = iter + 1;
    time = time + dt;
    timeVec(iter) = time;
    
    rms(iter) = 1/length(x) * sqrt(sum(x.^2));
    
    figure(1)
    subplot(2,1,1)
    plot(x,'.')
    ylim([-2 2])
    grid on
    title(sprintf('iter %d, time = %f', iter, time))
    
    subplot(2,1,2)
    plot(timeVec, rms)
    
    drawnow
    pause(0.01)
    

    E_h = Q*x;

    pos = primal.edgePos;
    vec = E_h / primal.edgeLength * primal.edgeDir;
    figure(2)
    quiver3(pos(:,1),  pos(:,2),  pos(:,3), ...
        vec(:,1), vec(:,2), vec(:,3),0)


    
    if iter > 10
        delta_rms = max(abs(rms(end-10 : end)) - mean(rms(end-10:end)));
        if delta_rms  < 1e-5
            break
        end
    end
end

%%

E_h = Q*x;

pos = primal.edgePos;
vec = E_h / primal.edgeLength * primal.edgeDir;
figure(2)
quiver3(pos(:,1),  pos(:,2),  pos(:,3), ...
    vec(:,1), vec(:,2), vec(:,3),0)


