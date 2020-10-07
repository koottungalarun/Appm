% function test_DEC()

% Solve the Gauss equation div(D) = rho with electrostatic potential
% boundary conditions, with d.o.f. given by the potential at vertices.

voltage_tA = 1;
voltage_tB = 2;

primal = readMesh('primal');
dual = readMesh('dual');

G = primal.e2v;
Sdual = dual.c2f;
Npf = size(primal.f2e,1);
Npe = size(G,1);
Npv = size(G,2);
values = dual.faceArea(1:Npe) ./ primal.edgeLength;
Meps = spdiags(values, 0, Npe, Npe);

S = Sdual(:,1:Npe);
assert(nnz(-G' - S) == 0) % assert that -G' == S


%% Boundary and Terminal Vertices
primal.isTerminalVertex = primal.vertexType == 2;
idx = true(size(primal.isTerminalVertex));
plot3(primal.vc(idx,1), primal.vc(idx,2), primal.vc(idx,3), '.')
hold on
idx = primal.isTerminalVertex;
plot3(primal.vc(idx,1), primal.vc(idx,2), primal.vc(idx,3), 'go')
idx = primal.isBoundaryVertex;
plot3(primal.vc(idx,1), primal.vc(idx,2), primal.vc(idx,3), 'r.')
hold off

%% Set equation and boundary conditions
A = S * Meps * (-G);

Npbv = sum(primal.isBoundaryVertex);
Nptv = sum(primal.isTerminalVertex);
phiA = voltage_tA * ones(Nptv/2,1);
phiB = voltage_tB * ones(Nptv/2,1);

phi = zeros(Npv,1);
idx = primal.isTerminalVertex;
phi(idx) = [phiA; phiB];

plot(primal.vc(:,3), phi, '*')

plot(phi,'*')
shg

spy(A)
hold on
line([1 size(A,2)], Nptv * [1 1], 'Color', 'r')
line(Nptv * [1 1], [1 size(A,1)], 'Color', 'r')
hold off


%% Solve system of equations
idx_d = false(size(phi));
idx_d(primal.isTerminalVertex) = true;
idx_f = ~idx_d;

A_free = A;
A_free(~idx_f, :) = [];
A_free(:,~idx_f) = [];
A_D = A(:, idx_d);

spy(A, 'g.')
hold on
spy(A_free, 'r.')
hold off
axis auto

rhs_D = A_D * phi(idx_d);
rhs = zeros(size(phi));
rhs_bc = rhs - rhs_D;

xfree = A_free \ rhs_bc(idx_f);

x = [xfree;
    phi(idx_d)];

%%
plot(primal.vc(:,3), x, '.')

E = -G * x;
plot(E, '.')

pos = primal.edgePos;
vec = E ./ primal.edgeLength .* primal.edgeDir;
quiver3(pos(:,1), pos(:,2), pos(:,3), ...
    vec(:,1), vec(:,2), vec(:,3), 0)
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

% end