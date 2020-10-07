clear
close all
clc

primalMesh = readMesh('primal');
dualMesh = readMesh('dual');

nPv = size(primalMesh.e2v, 2); % number of primal vertices
nPvi = sum(primalMesh.vertexType == 3); % inner vertices
nPvb = sum(primalMesh.vertexType <= 2); % boundary vertices
nPvt = sum(primalMesh.vertexType == 2); % terminal vertices

nPe = size(primalMesh.e2v, 1); % Number of primal edges
nPei = sum(primalMesh.edgeType <= 1); % number of primal edges in interior
nPeb = sum(primalMesh.edgeType == 2); % number of primal edges on boundary

nPf = size(primalMesh.f2e,1);
nPfi = sum(primalMesh.isFaceBoundary == 0);

Mnu =  spdiags(dualMesh.edgeLength(1:nPf) ./ primalMesh.faceArea, 0, nPf, nPf);
Mnu00 = Mnu(1:nPfi, 1:nPfi);

Meps = spdiags(dualMesh.faceArea(1:nPe) ./ primalMesh.edgeLength, 0, nPe, nPe);
G = primalMesh.e2v;
X = vertcat(sparse(nPvi, nPvb), speye(nPvb));
GX = -G * X;
Q = horzcat(vertcat(speye(nPei), sparse(nPeb, nPei)), GX);
C = primalMesh.f2e;
C00 = C(1:nPfi, 1:nPei);
T = [C00 sparse(nPfi, nPvb)];

M1 = Q' * Meps * Q;
M2 = T' * Mnu00 * T;
dt = 1;
M = M1 + dt^2 * M2;
spy(M)

x = zeros(nPei + nPvb,1);
assert(length(x) == size(M,1));
idx_d = false(size(x));
idx_d(end - ( (1:nPvt) - 1)) = true;
idx_f = ~idx_d;

M_f = M(idx_f, idx_f);
M_d = M(:, idx_d);

condest(M)
condest(M_f)

x_d = zeros(nPvt,1);
x_d(1:end/2) = 1;
x_d((end/2 + 1) : end) = 2;

rhs = zeros(size(x));
rhs = rhs - M_d * x_d;
rhs_f = rhs(idx_f);

x_f = M_f \ rhs_f;

figure(1)
x = [x_f; x_d];
plot(x, '.')

E_h = Q*x;
pos = primalMesh.edgePos;
vec = primalMesh.edgeDir .* repmat(E_h ./ primalMesh.edgeLength, 1, 3);
figure(2)
quiver3(pos(:,1), pos(:,2), pos(:,3), ...
    vec(:,1), vec(:,2), vec(:,3), 0)
axis equal
return
%%


% Q = readSparseMatrix('Q.dat');
C = readSparseMatrix('C.dat');
Mnu = readSparseMatrix('Mnu.dat');
Meps = readSparseMatrix('Meps.dat');
G = primalMesh.e2v;

% Number of vertices on boundary 
% (terminal vertices + free surface vertices)
nPvb = sum(primalMesh.vertexType == 1) + sum(primalMesh.vertexType == 2);

% Number of vertices in interior
nPvi = sum(primalMesh.vertexType == 3);

% Extension of interior and boundary vertices to all vertices
X = vertcat( sparse(nPvi, nPvb), speye(nPvb));

% Number of interior edges + interior-to-boundary edges
nPei = sum(primalMesh.edgeType == 0) + sum(primalMesh.edgeType == 1);

% Number of primal edges
nPe = size(G,1);

Q = horzcat( vertcat(speye(nPei), sparse(nPe - nPei, nPei)), G*X);

spy(Q)

condest(Q' * Q)

return


isTerminalVertex = primalMesh.vertexType == 2;

vc = primalMesh.vc;
idx = isTerminalVertex;
plot3(vc(idx,1), vc(idx,2), vc(idx,3), 'r.')
hold on
idx = primalMesh.vertexType == 1;
plot3(vc(idx,1), vc(idx,2), vc(idx,3), 'g.')
idx = primalMesh.vertexType == 3;
plot3(vc(idx,1), vc(idx,2), vc(idx,3), 'b.')
hold off

%%
% 
% phi = vc(:,3);
% % plot(phi)
% 
% E_h = G * phi;
% % plot(E_h)
% 
% temp = E_h ./ primalMesh.edgeLength;
% % plot(temp);
% 
% dir = primalMesh.edgeDir;
% dir = dir .* repmat(temp, 1, 3);
% pos = primalMesh.edgePos;
% % quiver3(pos(:,1), pos(:,2), pos(:,3), dir(:,1), dir(:,2), dir(:,3), 0)
% % axis equal
%%

nPei = sum(primalMesh.edgeType == 0);
nPeib = sum(primalMesh.edgeType == 1);
nPvb = sum(primalMesh.vertexType == 1);
nPvt = sum(primalMesh.vertexType == 2);

N = sum([nPei nPeib nPvb nPvt]);
assert(size(Q,1) == length(primalMesh.edgeLength))
assert(N == size(Q,2));


x = zeros(N,1);
idx_f = (1:N) <= (nPei + nPeib + nPvb);
idx_d = ~idx_f;

% electric potetial on terminal vertices
phi_t = [...
    1 * ones(nPvt/2,1);
    2 * ones(nPvt/2,1)
    ];
x_d = phi_t;

assert(sum(idx_d) == length(phi_t))

A = Q' * Meps * Q;

rhs = zeros(N,1);
A_d = A(:,idx_d);
rhs = rhs - A_d * x_d;
rhs_f = rhs(idx_f);
A_f = A(idx_f, idx_f);
x_f = A_f \ rhs_f;


