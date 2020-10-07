% Solve div(D) = rho with electrostatic potential boundary conditions and 
% d.o.f. given by inner edges and boundary vertex potentials.

clear
clc

primal = readMesh('primal');
dual = readMesh('dual');

Npv = size(primal.e2v,2);
Npe = size(primal.e2v,1);
Npf = size(primal.f2e,1);
G = primal.e2v;
Meps = spdiags(dual.faceArea(1:Npe) ./ primal.edgeLength, 0, Npe, Npe);
isTerminalVertex = primal.vertexType == 2;

%% Plot edges
pos = primal.edgePos;
vec = primal.edgeDir;

idx = primal.edgeType == 0; % inner edges
quiver3(pos(idx,1), pos(idx,2), pos(idx,3), ...
    vec(idx,1), vec(idx,2), vec(idx,3), 0, 'r')
hold on
idx = primal.edgeType == 1; % inner-to-boundary edges
quiver3(pos(idx,1), pos(idx,2), pos(idx,3), ...
    vec(idx,1), vec(idx,2), vec(idx,3), 0, 'g')
idx = primal.edgeType == 2;% boundary edges
quiver3(pos(idx,1), pos(idx,2), pos(idx,3), ...
    vec(idx,1), vec(idx,2), vec(idx,3), 0, 'b')
hold off

%%
idx = primal.vertexType == 1;
pos = primal.vc;
plot3(pos(idx,1), pos(idx,2), pos(idx,3), '.r')
hold on
idx = primal.vertexType == 2;
plot3(pos(idx,1), pos(idx,2), pos(idx,3), '.g')
idx = primal.vertexType == 3;
plot3(pos(idx,1), pos(idx,2), pos(idx,3), '.b')
hold off

%%
Npvb = sum(primal.vertexType == 1 | primal.vertexType == 2);
Npvi = sum(primal.vertexType == 3);
assert((Npvb + Npvi) == Npv)
Npei = sum(primal.edgeType <= 1);


Q = sparse(Npe, Npei + Npvb);
for k = 1 : Npe
    if primal.edgeType(k) <= 1
        Q(k,k) = 1;
    end
%     continue
    if primal.edgeType(k) >= 1
        [~,col,val] = find(primal.e2v(k,:));
        for m = 1 : length(col)
            vIdx = col(m);
            if primal.vertexType(vIdx) <= 2
                p = Npei + vIdx - Npvi;
                assert(p > 0)
                assert(p <= size(Q,2))
                Q(k, p) = -1 * val(m);
            end
        end
    end
end

spy(Q)
shg

%%
G = primal.e2v;
phi = primal.vc(:,3);
E = -G * phi;

vec = E ./ primal.edgeLength .* primal.edgeDir;
pos = primal.edgeCenter;
quiver3(pos(:,1), pos(:,2), pos(:,3), ...
    vec(:,1), vec(:,2), vec(:,3), 0)

%%
Sdual = dual.c2f;
S = Sdual(:,1:Npe);
D = Meps * E;
% plot(E ./ primal.edgeLength .* dual.faceArea(1:Npe),'.')
plot(D,'.')
pos = dual.faceCenter(1:Npe,:);
vec = D ./ dual.faceArea(1:Npe) .* dual.faceNormal(1:Npe,:);
quiver3(pos(:,1), pos(:,2), pos(:,3), ...
    vec(:,1), vec(:,2), vec(:,3))

Sdual = dual.c2f;
S = Sdual(:,1:Npe);
rho = S * D;
idx = abs(rho) > 1e-10;
plot3(primal.vc(idx,1), primal.vc(idx,2), primal.vc(idx,3), 'o')

%%
E_i = E(1:Npei);
phi_b = phi(primal.isBoundaryVertex);
x = [E_i; phi_b];
E2 = Q*x;

plot(E2 ./ primal.edgeLength, '.')
vec = E2 ./ primal.edgeLength .* primal.edgeDir;
pos = primal.edgePos;
quiver3(pos(:,1), pos(:,2), pos(:,3), ...
    vec(:,1), vec(:,2), vec(:,3))
axis equal

%%
% Sdual = dual.c2f;
% S = Sdual(:, 1:Npe);
A = Q' * Meps * Q;
spy(A)


%%
% boundary vertices: vertexType=1 (boundary) or 2 (terminal)
% inner edges: edgeType=0 or edgeType=1
Npvt = sum(primal.vertexType == 2);
assert(mod(Npvt,2) == 0);
phi1 = 1 * ones(Npvt/2,1);
phi2 = 2 * ones(Npvt/2,1);
phi_t = [phi1; phi2];
phi_b = [
    zeros(Npvb - Npvt, 1);
    phi_t];
assert(length(phi_b) == Npvb)
E_i = zeros(Npei,1);
x = [E_i; phi_b];
assert(length(x) == size(A,2))
plot(x, '.')

%% Solve system of equations
% 
idx_d = false(size(x));
idx_d(end + 1 - (1:Npvt)) = true;
idx_f = ~idx_d;

Afree = A(idx_f,idx_f);
A_D = A(:, idx_d);
rhs = zeros(size(A,1),1);
rhsfree = rhs - A_D * phi_t; %x(idx_d);
xfree = Afree \ rhsfree(idx_f);

x(idx_f) = xfree;
plot(x,'.')


E = Q*x;
plot(E,'.')

% vec = E ./ primal.edgeLength .* primal.edgeDir;
% pos = primal.edgeCenter;
% quiver3(pos(:,1), pos(:,2), pos(:,3), ...
%     vec(:,1), vec(:,2), vec(:,3), 0)


