clear
clc

primal = readMesh('primal');
dual = readMesh('dual');


%% Number of mesh items (vertices, edges, faces)
Npv = size(primal.e2v,2); % no. of primal vertices
Npvi = sum(primal.vertexType == 3); % no. of primal inner vertices
Npvb = sum(primal.vertexType <= 2); % no. of primal boundary vertices
Npvt = sum(primal.vertexType == 2); % no. of primal terminal vertices

Npe = size(primal.e2v,1); % no. of primal edges
Npei = sum(primal.edgeType <= 1); %inner edges
Npeb = sum(primal.edgeType == 2); % boundary edges

Npf = size(primal.f2e,1); % no. of primal faces
Npfi = sum(primal.isFaceBoundary == 0);
Npfb = sum(primal.isFaceBoundary == 1);

%% Define faces for current
pos = dual.faceCenter(1:Npe,:);
fn  = dual.faceNormal(1:Npe,:);

idx = false(Npe,1);
idx(1:Npei) = true;

zUnitVec = repmat([0 0 1], Npe, 1);
idx2 = vecnorm(cross(fn', zUnitVec'))' < 0.1;
pos_2d = pos(1:Npe,1:2);
idx3 = vecnorm(pos_2d')' < 0.35;
idx = idx & idx2 & idx3;

currentIdx = idx;

vec = fn;
quiver3(pos(idx,1), pos(idx,2), pos(idx,3), ...
    vec(idx,1), vec(idx,2), vec(idx,3))


%% Plot vertex mesh with vertex types
pos = primal.vc;

idx = primal.vertexType == 1; % boundary, not terminal
plot3(pos(idx,1), pos(idx,2), pos(idx,3), 'r.')
hold on
idx = primal.vertexType == 2; % boundary, terminal
plot3(pos(idx,1), pos(idx,2), pos(idx,3), 'g.')
idx = primal.vertexType == 3; % inner
plot3(pos(idx,1), pos(idx,2), pos(idx,3), 'b.')
hold off


%% Plot edges and edge type

pos = primal.edgePos;
vec = primal.edgeDir;

idx = primal.edgeType == 0;
quiver3(pos(idx,1), pos(idx,2), pos(idx,3), ...
    vec(idx,1), vec(idx,2), vec(idx,3), 0, 'Color', 'r')
hold on
idx = primal.edgeType == 1;
quiver3(pos(idx,1), pos(idx,2), pos(idx,3), ...
    vec(idx,1), vec(idx,2), vec(idx,3), 0, 'Color', 'g')
idx = primal.edgeType == 2;
quiver3(pos(idx,1), pos(idx,2), pos(idx,3), ...
    vec(idx,1), vec(idx,2), vec(idx,3), 0, 'Color', 'b')
hold off


%% Voltage boundary condition on terminal vertices
% phi1 = 1;
% phi2 = 2;

% phi_b = [zeros(Npvb - Npvt,1);
%     phi_t];


% E_i = zeros(Npei,1);
% x = [E_i; phi_b];
%% Operator C

C = primal.f2e;
Cdual = dual.f2e;

C00     = C(1:Npfi, 1:Npei);
C00dual = Cdual(1:Npei, 1:Npfi);

figure(1)
clf
spy(C)
hold on
spy(C00, 'r')
hold off
axis tight

figure(1)
subplot(1,2,1)
spy(C00')
subplot(1,2,2)
spy(C00dual)
clf

T = C00' - C00dual;
assert( nnz(T) == 0 )
clear T

%% Operator G, X, and Q
G = primal.e2v;
X = vertcat(sparse(Npvi, Npvb), speye(Npvb, Npvb));
GX = G*X;
spy(GX)

% T = [speye(Npei);
%     sparse(Npeb, Npei)];
T = vertcat(speye(Npei), sparse(Npeb, Npei));
Q = horzcat(T, GX);
% clear T

spy(Q)

%% Operator Meps and Mnu
values = dual.faceArea(1:Npe) ./ primal.edgeLength;
Meps = spdiags(values, 0, Npe, Npe);

values = dual.edgeLength(1:Npf) ./ primal.faceArea;
Mnu = spdiags(values, 0, Npf, Npf);

Mnu00 = Mnu(1:Npfi, 1:Npfi);

T = [C00 sparse(Npfi, Npvb)];

%% Simulation parameters
iterMax = 15000;
timeMax = 100;
dt = 1e-0;
plotIter = min([1, floor(1/dt), iterMax / 20]);
iter = 1;
t = 0;

phi1    = @(tt) 0;
phi2    = @(tt) 0; %1 - exp(-tt / 4);
current = @(tt) 10 * 1/2 * (1 + tanh( (tt - 50) / 10 ));

J = zeros(Npe,1);
%% System of equations
A = Q'*Meps*Q;
B = T' * Mnu00 * T;
M = A + dt^2 * B;

x = zeros(Npei + Npvb,1);
assert(length(x) == size(M,2))

idx_d = false(size(x));
idx_d( end + 1 - (1:Npvt) ) = true;
idx_f = ~idx_d;

M_D = M(:,idx_d);
M_f = M(idx_f, idx_f);

%%
rhs = zeros(size(x));
x_mm1 = zeros(size(x));
x_m = zeros(size(x));
current_dt = 0;
while t <= timeMax && iter <= iterMax
    % Solve for state vector x at timestep x+1
    rhs = zeros(size(x));
    
    currentVec(iter) = current(t);
    if iter > 1
        current_dt(iter) = (currentVec(iter) - currentVec(iter-1)) / dt;
        rhs(currentIdx) = current_dt(iter);
    end
    rhs = A * (2*x_m - x_mm1) + dt^2 * rhs;
    
    phi_t = [
        phi1(t) * ones(Npvt/2,1); 
        phi2(t) * ones(Npvt/2,1)
    ];

    x_d = phi_t;
    rhs = rhs - M_D * x_d;
    rhs_f = rhs(idx_f);
    x_f = M_f \ rhs_f;
    x = [x_f; x_d];
    E_h = Q*x;
    
    % update of state vectors at timestep m and m-1 (m minus 1)
    x_mm1 = x_m;
    x_m = x;
    t = t + dt;
    iter = iter + 1;

    rms(iter) = sqrt(sum( (x_m - x_mm1).^2 )) / length(x_m);
    emean(iter) = mean(E_h);
    time(iter) = t;
    
    % Plot data
    if mod(iter, plotIter) == 0
        figure(1)
        subplot(3,2,[1 3 5])
        pos = primal.edgeCenter;
        vec = E_h ./ primal.edgeLength .* primal.edgeDir;
        figure(1)
        quiver3(pos(:,1), pos(:,2), pos(:,3), ...
            vec(:,1), vec(:,2), vec(:,3), 0)
%         zlim([-0.2 1.2])
        %     axis equal
        title(iter)
        drawnow
        
        subplot(3,2,2)
        semilogy(time, rms)
        title('rms')
        grid on
        ylim([1e-16 1])
        
        subplot(3,2,4)
        plot(time, emean)
        title('E mean')
        grid on
        
        subplot(3,2,6)
        if ~isempty(current_dt)
            plot(current_dt)
        end
        title('dJ/dt')
        grid on
    end
end


%%
