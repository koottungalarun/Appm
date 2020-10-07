clear
clc

primal = readMesh('primal');
dual = readMesh('dual');

G = primal.e2v;
C = primal.f2e;

nPv = size(G,2);
nPvb = sum(primal.vertexType <= 2);
nPvi = sum(primal.vertexType == 3);
nPvt = sum(primal.vertexType == 2);

nPe = size(G,1);
nPei = sum(primal.edgeType <= 1);
nPeb = sum(primal.edgeType == 2);

nPf = size(C,1);
nPfi = sum(primal.isFaceBoundary == 0);
nPfb = sum(primal.isFaceBoundary == 1);


%%
x_h = zeros(nPei + nPvb,1);
X = vertcat(sparse(nPvi, nPvb), speye(nPvb));
GX = G * X;
Q = horzcat(vertcat(speye(nPei), sparse(nPeb, nPei)), -GX);

C_00 = C(1:nPfi, 1:nPei);
P = horzcat(C_00, sparse(nPfi, nPvb));

Meps = spdiags(dual.faceArea(1:nPe) ./ primal.edgeLength, 0, nPe, nPe);
Mmu = spdiags(primal.faceArea(1:nPfi) ./ dual.edgeLength(1:nPfi), 0, nPfi, nPfi);

tau = 0.1;
A = Q' * Meps * Q + tau^2 * P' * Mmu * P;
idx_d = false(size(A,2),1);
idx_d(end + 1 - (1:nPvt/2)) = true;
idx_f = ~idx_d;
A_F = A(idx_f,idx_f);
A_D = A(:,idx_d);
condest(A_F)

phi0 = @(t) 0;
phi1 = @(t) (t > 0) .* (1 * 0.5 * (1 + tanh( (t - 5) / 1)));
phiT = @(t) [phi0(t) * ones(nPvt/2,1); phi1(t) * ones(nPvt/2,1); ];

iter = 1;
time = tau;
timeMax = 50;
iterMax = min(50000,ceil(timeMax/tau));

idxE = primal.edgeCenter(:,3) == 0.5 & vecnorm(primal.edgeCenter(:,1:2)')' == 0 & vecnorm(primal.edgeDir(:,1:2)')' < 0.1;
assert(sum(idxE) == 1)

x_prev = x_h;
temp = zeros(nPe,1);

% J_src = @(t) exp(-(t - 5).^2 / 0.5) * (ones(nPe,1) ./ dual.faceArea(1:nPe) .* primal.edgeLength(1:nPe));
J_src = @(t) exp(-(t - 10).^2 / 0.5) * 1e-6 * [(ones(nPei,1) ./ dual.faceArea(1:nPei) .* primal.edgeLength(1:nPei)); zeros(nPe-nPei,1)];

while iter <= iterMax && time <= timeMax
    timeVec(iter) = time;
    
    x_prev2 = x_prev;
    x_prev = x_h;
    
%     x_d = phiT(time);
    x_d = phi0(time) * ones(nPvt/2,1);
    rhs = zeros(size(x_h));
    rhs = rhs + Q'*Meps*Q * (2*x_prev - x_prev2);
    rhs = rhs - A_D * x_d;
    rhs = rhs - Q' * (J_src(time + tau) - J_src(time));
    
    x_f = A_F \ rhs(idx_f);
    x_h = [x_f; x_d];
    
    E_h = Q * x_h;
    rms(iter) = 1/length(E_h) * sqrt(sum(E_h(:).^2));
    dataE(iter) = E_h(idxE) / primal.edgeLength(idxE);
    
    if mod(iter,10) == 0
    pos = primal.edgeCenter;
    vec = E_h ./ primal.edgeLength .* primal.edgeDir;
    figure(1)
    subplot(2,2,[1 3])
    quiver3(pos(:,1), pos(:,2), pos(:,3), vec(:,1), vec(:,2), vec(:,3))
    zlim([-0.2 1.2])
    subplot(2,2,2)
    plot(timeVec,rms)
    grid on
    subplot(2,2,4)
    plot(timeVec,dataE)
    grid on
    title(iter)
    drawnow
    end
    
    time = time + tau;
    iter = iter + 1;    
end



