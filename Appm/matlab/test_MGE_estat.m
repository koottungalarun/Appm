% Test Maxwell Grid Equations with electrostatic boundary conditions

clear
clc

%% Read mesh data
primal = readMesh('primal');
dual = readMesh('dual');

%% Define primal boundary vertex type
primal.isTerminalVertex = primal.isBoundaryVertex & vecnorm(primal.vc(:,1:2)')' < 1.5 ...
    & (primal.vc(:,3) == 0 | primal.vc(:,3) == 1);
boundaryVertexType = zeros(size(primal.isBoundaryVertex));

idx = ~(primal.isBoundaryVertex);
boundaryVertexType(idx) = 0;
idx = primal.isBoundaryVertex & ~(primal.isTerminalVertex);
boundaryVertexType(idx) = 1;
idx = primal.isTerminalVertex;
boundaryVertexType(idx) = 2;

primal.boundaryVertexType = boundaryVertexType;
clear boundaryVertexType

% Plot primal vertices
idx = primal.boundaryVertexType == 0;
figure(1)
plot3(primal.vc(idx,1), primal.vc(idx,2), primal.vc(idx,3), 'b.')
hold on
idx = primal.boundaryVertexType == 1;
plot3(primal.vc(idx,1), primal.vc(idx,2), primal.vc(idx,3), 'r.')
idx = primal.boundaryVertexType == 2;
plot3(primal.vc(idx,1), primal.vc(idx,2), primal.vc(idx,3), 'g.')
hold off
xlabel('x')
ylabel('y')
zlabel('z')
axis equal

%% Plot primal edges
plot(primal.edgeType)

idx = primal.edgeType == 0;
quiver3(primal.edgePos(idx,1), primal.edgePos(idx,2), primal.edgePos(idx,3), ...
    primal.edgeDir(idx,1), primal.edgeDir(idx,2), primal.edgeDir(idx,3), 0, 'Color', 'b')
hold on
idx = primal.edgeType == 1;
quiver3(primal.edgePos(idx,1), primal.edgePos(idx,2), primal.edgePos(idx,3), ...
    primal.edgeDir(idx,1), primal.edgeDir(idx,2), primal.edgeDir(idx,3), 0, 'Color', 'r')
idx = primal.edgeType == 2;
quiver3(primal.edgePos(idx,1), primal.edgePos(idx,2), primal.edgePos(idx,3), ...
    primal.edgeDir(idx,1), primal.edgeDir(idx,2), primal.edgeDir(idx,3), 0, 'Color', 'g')
hold off

%% Define extension mapping E = X * x, with x = [E_i, phi_b]
Npv = size(primal.vc,1); % number of primal vertices
Npbv = sum(primal.isBoundaryVertex); % number of primal boundary vertices
Npiv = Npv - Npbv; % number of primal inner vertices
X = [sparse(Npiv, Npbv);
    speye(Npbv)];

% Define electric potential at boundary vertices
idx = primal.isBoundaryVertex;
phiB = primal.vc(idx,2); % potential = z-coordinate

Npe = size(primal.e2v,1); % number of primal edges
Npbe = sum(primal.edgeType >= 1); % number of primal boundary edges and those connecting to interior
Npie = Npe - Npbe; % number of primal inner edges
Npe0 = sum(primal.edgeType == 0); % number of primal edges of type 0 (inner)
assert(Npie == Npe0)
Npe1 = sum(primal.edgeType == 1); % number of primal edges of type 1 (inner-to-boundary)
Npe2 = sum(primal.edgeType == 2); % number of primal edges of type 2 (boundary)


% Gradient on primal mesh
G = primal.e2v;
% spy(G)
% hold on
% line((Npiv + 0.5) * [1 1], [0 Npe], 'Color', 'r')
% line([0 Npv], (Npie + 0.5) * [1 1], 'Color', 'r')
% hold off
% shg

phi = primal.vc(:,3);
E = -G * phi;
vec = E ./ primal.edgeLength .* primal.edgeDir;
quiver3(primal.edgePos(:,1), primal.edgePos(:,2), primal.edgePos(:,3), ...
    vec(:,1), vec(:,2), vec(:,3), 0)

%%
Q = sparse(Npe, Npe0 + Npe1 + Npbv);
for k = 1 : Npe
    if primal.edgeType(k) <= 1
        Q(k,k) = 1;
    end
    if primal.edgeType(k) >= 1
        [~,vIdx,incid] = find(primal.e2v(k,:));
        nv = length(vIdx);
        for m = 1 : nv
            if primal.isBoundaryVertex(vIdx(m))
                vOffset = vIdx(m) - Npiv;
                assert(vOffset > 0)
                p = Npe0 + Npe1 + vOffset;
                Q(k,p) = -incid(m);
            end
        end
    end
end
spy(Q)


%% Inclusion-/gradient operator
Ei = zeros(Npe0 + Npe1, 1);
% inclusion operator for inner edges
Incl_E = [speye(Npe0 + Npe1); sparse(Npe2, Npe0 + Npe1)]; 
E = Incl_E * Ei + G * X * phiB;


% inclusion-gradient operator for vector x (vector of d.o.f.)
% E = Q*x
Q_e = [Incl_E, -G*X];
assert(nnz(Q - Q_e) == 0)

x = [Ei; phiB];
assert(size(Q_e,2) == size(x,1))

%% Solve electrostatic equation div(D) = 0
values = dual.faceArea(1:Npe) ./ primal.edgeLength;
Meps = spdiags(values, 0, Npe, Npe);
A = Q' * Meps * Q;
% spy(A)

E = -G * phi;

idx_d = false(size(x));
idx_d((Npe0 + Npe1 + 1) : end) = true;
idx_f = ~idx_d;

A_D = A(:,idx_d);
rhs = zeros(size(x));
rhs_D = rhs - A_D * x(idx_d);

A_f = A(idx_f, idx_f);
x_f = A_f \ rhs_D(idx_f);


return

%%
idx = primal.edgeType <= 1;
assert(sum(idx) == Npe0 + Npe1)
E_i = E(idx);
phiB = phi(primal.isBoundaryVertex);
x = [E_i; phiB];
E = Q_e * x;
vec = E ./ primal.edgeLength .* primal.edgeDir;
quiver3(primal.edgePos(:,1), primal.edgePos(:,2), primal.edgePos(:,3), ...
    vec(:,1), vec(:,2), vec(:,3), 0)



%% Material laws
Npf = size(primal.f2e,1); % number of primal faces
Npbf = sum(primal.isFaceBoundary);
Npif = sum(~(primal.isFaceBoundary));
assert(Npbf + Npif == Npf)

values = dual.faceArea(1:Npe) ./ primal.edgeLength;
Meps = spdiags(values, 0, Npe, Npe);

values = primal.faceArea ./ dual.edgeLength(1:Npf);
Mnu = spdiags(values, 0, Npf, Npf);

%% Curl-operator for inner primal edges / inner primal faces
C = primal.f2e; % complete curl-operator on primal mesh
C_i = C(1:Npif, 1:Npie); % curl-operator on interor of primal mesh

Q_nu = [C_i sparse(Npif, Npe1 + Npbv)];
assert(size(Q_nu, 2) == size(x,1))


%% Solve reformulated ampere equation for inner faces
A =  Q_e' * Meps * Q_e;
B = Q_nu' * Mnu(1:Npif,1:Npif) * Q_nu;
assert(all(size(A) == size(B)))

dt = 1;
M = dt^2 * A + B;

x = zeros(Npe0 + Npe1 + Npbv, 1);
idx = false(size(x));

Nptv = sum(primal.isTerminalVertex); % number of primal boundary vertices on terminals
idx_d = false(size(x)); % index of fixed values (Dirichlet b.c.)
idx_f = false(size(x)); % index of free values 

idx_d(end + 1 - (1:Nptv)) = true; 
idx_f = ~idx_d;

assert(mod(Nptv,2) == 0)
phiA = 1 * ones(Nptv/2,1);
phiB = 2 * ones(Nptv/2,1);
assert(sum(idx_d) == length(phiA) + length(phiB))
x(idx_d) = [phiA; phiB];

rhs = zeros(size(x));

rhs_D = M(:, idx_d) * x(idx_d);
xfree = M(:, idx_f) \ (rhs - rhs_D);

x(idx_f) = xfree;
plot(x,'.')
% hold on
% plot(xfree, '.')
% hold off
shg

% Goal:
% (d/dt)^2 ( Q' M_eps Q x ) + (C' M_nu C) x = rhs
% ->
% (dt^2 * Q' Meps Q + C' Mnu * C) * x^(m+1) = dt^2 * (Q' * Meps * Q) *
% (2*x^m - x^(m-1)) + rhs

%%
E = Q_e * x;
evec = E ./ primal.edgeLength .* primal.edgeDir;
quiver3(primal.edgeCenter(:,1), primal.edgeCenter(:,2), primal.edgeCenter(:,3), ...
    evec(:,1), evec(:,2), evec(:,3), 0)