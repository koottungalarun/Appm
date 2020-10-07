clear
clc

%% Load data
tau = 1;
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


%% Setup operators
x_h = zeros(nPei + nPvb,1);
X = vertcat(sparse(nPvi, nPvb), speye(nPvb));
GX = G * X;
Q = horzcat(vertcat(speye(nPei), sparse(nPeb, nPei)), -GX);

C_00 = C(1:nPfi, 1:nPei);
P = horzcat(C_00, sparse(nPfi, nPvb));

Meps = spdiags(dual.faceArea(1:nPe) ./ primal.edgeLength, 0, nPe, nPe);
Mmu = spdiags(primal.faceArea(1:nPfi) ./ dual.edgeLength(1:nPfi), 0, nPfi, nPfi);

H_h = zeros(nPfi,1);
u_h = [H_h; x_h];

%% Setup system of equations for Crank-Nicholson scheme
A_11 = Mmu;
A_12 = 0.5 * tau * P;
A_21 = -0.5 * tau * P';
A_22 = Q' * Meps * Q;

B_11 = Mmu;
B_12 = -0.5 * tau * P;
B_21 = 0.5 * tau * P';
B_22 = Q' * Meps * Q;

A = [A_11, A_12;
    A_21, A_22];
B = [B_11, B_12;
    B_21, B_22];

assert(size(A,2) == size(u_h,1))
assert(size(B,2) == size(u_h,1))

%% Setup of current source term
Jz_src = dual.faceArea ...
    .* dot(dual.faceNormal', repmat([0 0 1], size(dual.faceNormal,1), 1)')' ...
    .* (vecnorm(dual.faceCenter(:,1:2)')' < 0.35);
J_src = zeros(nPe,1);
J_src(1:nPei) = Jz_src(1:nPei);
assert(size(J_src,1) == nPe);

% pos = dual.faceCenter;
% vec = dual.faceNormal;
% figure(1)
% clf
% quiver3(pos(:,1), pos(:,2), pos(:,3), vec(:,1), vec(:,2), vec(:,3))
% hold on
% idx = 1:nPe;
% vec = repmat(J_src ./ dual.faceArea(idx), 1, 3) .* dual.faceNormal(idx,:);
% quiver3(pos(idx,1), pos(idx,2), pos(idx,3), ...
%     vec(idx,1), vec(idx,2), vec(idx,3), 'LineWidth', 2)
% hold off

% rhs = [
%     zeros(nPfi,1);
%     -0.5 * tau * Q2' * J_src];

rhs = [zeros(nPfi,1);
    J_src(1:nPei);
    zeros(nPvb,1);
    ];
size(rhs,1) == size(A,1);

% J_edge = rhs(nPfi + 1 : end);
% plot(J_edge,'.')
% hold on
% line(nPei * [1 1], 2e-3 * [-1 1]', 'Color','k')
% hold off

%%
figure(1)
clf

idxD = false(size(u_h));
idxD(end+1 - (1:nPvt/2)) = true;
idxf = ~idxD;
Afree = A(idxf, idxf);
% condest(Afree)
AD = A(:,idxD);
u_d = zeros(sum(idxD),1);
u_h = zeros(size(u_h));

iterMax = 100;
visIter = 1;
for iter = 1 : iterMax
    
    u_h_previous = u_h;
    
    rhs = zeros(size(u_h));
    rhs = rhs + B * u_h_previous;
    rhs = rhs + [zeros(nPfi,1); -0.5 * tau * Q' * J_src];
%     rhs = rhs + [zeros(nPfi,1); -0.5 * tau * J_src(1:nPei); zeros(nPvb,1)];
    rhs = rhs - AD * u_d;
    rhsf = rhs(idxf);
    u_f = Afree \ rhsf;
    u_h = [u_f; u_d];
    
    H_h = u_h(1:nPfi);
    x_h = u_h(nPfi + 1 : end);
    E_h = Q*x_h;
    
    phiB = x_h(end + 1 - (1:nPvt));
    phiB_max(iter) = max(abs(phiB));
    E_h_max(iter) = max(abs(E_h));
    H_h_max(iter) = max(abs(H_h));
    
    if mod(iter,visIter) == 0 || iter == iterMax
        figure(1)
        pos = primal.edgeCenter;
        vec = repmat(E_h ./ primal.edgeLength, 1, 3) .* primal.edgeDir;
        subplot(1,2,1)
        quiver3(pos(:,1), pos(:,2), pos(:,3), vec(:,1), vec(:,2), vec(:,3), 1)
        grid on
        xlabel('x')
        ylabel('y')
        zlabel('z')
        
        
        subplot(1,2,2)
        plot(E_h,'.')
%         ylim(3*[-1 1])
        grid on
        ylabel('E_h')
        title(sprintf('iter %d', iter))
        drawnow

        figure(2)
        idx = 1 : nPfi;
        pos = dual.edgeCenter(idx,:);
        vec = repmat(H_h ./ dual.edgeLength(idx), 1, 3) .* dual.edgeDir(idx,:);
        subplot(1,2,1)
        quiver3(pos(:,1), pos(:,2), pos(:,3), vec(:,1), vec(:,2), vec(:,3))
        grid on
        xlabel('x')
        ylabel('y')
        zlabel('z')        

        subplot(1,2,2)
        plot(H_h,'.')
        ylim(0.1*[-1 1])
        grid on
        ylabel('H_h')
%         title(sprintf('iter %d', iter))
        drawnow
        
        figure(3)
        plot(phiB,'.')
        title('phiB')
        grid on
        drawnow
        
        figure(4)
        subplot(3,1,1)
        plot(phiB_max)
        grid on
        
        subplot(3,1,2)
        plot(E_h_max)
        grid on
        
        subplot(3,1,3)
        plot(H_h_max)
        grid on
        
        xlabel('iter')
        
    end
end
return
%%
Mnu = spdiags(dual.edgeLength(1:nPf) ./ primal.faceArea, 0, nPf, nPf);
T = Q'*C'*Mnu*C*Q;
spy(T)
% T2 = T(1:nPei, 1:nPei);
% spy(T2)
% condest(T2)

idxd = false(size(x_h));
idxd((end - nPvb) : end) = true;
idxf = ~idxd;
Tfree = T(idxf, idxf);
spy(Tfree)
condest(Tfree)