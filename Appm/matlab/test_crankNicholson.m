% Test for Crank-Nicholson time integration

clear
clc

%%
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
Ei = zeros(nPei,1);
phiB = zeros(nPvb,1);

zUnit = [0 0 1];
a = primal.edgeDir(1:nPei,:) ./ primal.edgeLength(1:nPei);
b = repmat(zUnit,nPei,1);
temp = dot(a',b')';
clear a b
plot(temp,'.')
idx = abs(temp) > 0.01;
Ei(idx) = primal.edgeLength(idx);

idx = primal.vertexType <= 2;
phiB = 1 - primal.vc(idx,3);

x_h = [Ei; phiB];
assert(all(size(x_h) == [nPei + nPvb,1]));

%%
X = vertcat(sparse(nPvi, nPvb), speye(nPvb));
GX = G * X;
Q = horzcat(vertcat(speye(nPei), sparse(nPeb, nPei)), -GX);
E_h = Q*x_h;

% Plot primal edge vectors to check if data layout is correct
pos = primal.edgeCenter;
vec = E_h ./ primal.edgeLength .* primal.edgeDir;
scale = 0;
% plot vertices
idx = primal.edgeType == 0;
plot3(primal.vc(:,1), primal.vc(:,2), primal.vc(:,3), '.')
% plot edges
hold on
quiver3(pos(idx,1), pos(idx,2), pos(idx,3), ...
    vec(idx,1), vec(idx,2), vec(idx,3), scale, 'Color', 'r')
idx = primal.edgeType == 1;
quiver3(pos(idx,1), pos(idx,2), pos(idx,3), ...
    vec(idx,1), vec(idx,2), vec(idx,3), scale, 'Color', 'g')
idx = primal.edgeType == 2;
quiver3(pos(idx,1), pos(idx,2), pos(idx,3), ...
    vec(idx,1), vec(idx,2), vec(idx,3), scale, 'Color', 'b')
hold off
clear idx

%%
C_00 = C(1:nPfi, 1:nPei);
P = horzcat(C_00, sparse(nPfi, nPvb));

%%
Meps = spdiags(dual.faceArea(1:nPe) ./ primal.edgeLength, 0, nPe, nPe);
Mmu = spdiags(primal.faceArea(1:nPfi) ./ dual.edgeLength(1:nPfi), 0, nPfi, nPfi);

pos = primal.edgePos;
vec = primal.edgeDir;
idx = true(nPe,1);
idx = vecnorm(pos(:,1:2)')' < 0.35;
quiver3(pos(idx,1), pos(idx,2), pos(idx,3), ...
    vec(idx,1), vec(idx,2), vec(idx,3), 0)

values = ones(nPe,1);
values(idx) = 10;
Msigma = spdiags(values, 0, nPe, nPe);

H_h = zeros(size(Mmu,2),1);


%%
tau = 0.1;
timeMax = 50;
iterMax = 5000;

T{1,1} = Mmu;
T{1,2} =  1/2 * tau * P;
T{2,1} = -1/2 * tau * P';
T{2,2} = Q'*Meps*Q + 1/2 * tau * Q' * Msigma * Q;

% If matrix sizes do not match, show the data structure
if (size(T{1,1},1) ~= size(T{1,2},1) || size(T{2,1},1) ~= size(T{2,2},2) || ...
        size(T{1,1},2) ~= size(T{2,1},2) || size(T{1,2},2) ~= size(T{2,2},2))
    disp(T)
    warning('Matrix sizes do not match')
end
A = [
    T{1,1}, T{1,2};
    T{2,1}, T{2,2}];

T{1,1} = Mmu;
T{1,2} = -1/2 * tau * P;
T{2,1} =  1/2 * tau * P';
T{2,2} = Q' * Meps * Q - 1/2 * tau * Q' * Msigma * Q;
B = [
    T{1,1}, T{1,2};
    T{2,1}, T{2,2}];

if (size(T{1,1},1) ~= size(T{1,2},1) || size(T{2,1},1) ~= size(T{2,2},2) || ...
        size(T{1,1},2) ~= size(T{2,1},2) || size(T{1,2},2) ~= size(T{2,2},2))
    disp(T)
    warning('Matrix sizes do not match')
end

u = [H_h; x_h];

% Check if data has correct size
assert(size(A,2) == size(u,1))
assert(size(B,2) == size(u,1))
assert(all(size(A) == size(B)))

%% Solve systems
iter = 1;
time = 0;
u = zeros(size(u));

idxd = false(size(u));
idxd(end + 1 - (1:nPvt)) = true;
idxf = ~idxd;

Af = A(idxf,idxf);
Ad = A(:,idxd);

phi0 = @(t) 0;
phi1 = @(t) (t > 0) * 10 * 0.5 * (1 + tanh( (t-5) ));
phiT = @(t) [
    phi0(t) * ones(nPvt/2, 1);
    phi1(t) * ones(nPvt/2, 1)];

ud = phiT(time);
assert(all(size(ud) == [nPvt,1]))
assert(sum(idxd) == nPvt)
assert(all(size(ud) == size(u(idxd))))

dt = tau;
timeVec = [];
while (time <= timeMax && iter <= iterMax)
    timeVec(end+1) = time;
    u_prev = u;
    ud = phiT(time);
    rhs = B * u_prev - Ad * ud;
    uf = Af \ rhs(idxf);
    u(idxf) = uf;
    u(idxd) = ud;
    
    rms(iter) = sqrt( sum( (u - u_prev).^2 ) ) / length(u);
    
    %     figure(1)
    %     subplot(2,1,1)
    %     plot(u,'.')
    %     grid on
    %     title(sprintf('Time = %d', time))
    %
    %     subplot(2,1,2)
    %     semilogy(timeVec,rms)
    %     title('rms')
    %     grid on
    %     drawnow
    
    
    assert(nPfi + nPei + nPvb == length(u))
    H_h = u(1:nPfi);
    H_h_prev = u_prev(1:nPfi);
    x_h = u(nPfi + (1 : (nPei + nPvb)));
    x_h_prev = u_prev(nPfi + (1 : (nPei + nPvb)));
    
    %     Gb = G(end+1-(1:nPeb), end+1-(1:nPvb));
    %     Q_bv = [
    %         speye(nPei), sparse(nPei, nPvb)
    %         sparse(nPeb, nPei), -Gb
    %         ];
    % x_h = zeros(nPei + nPvb,1);
    % x_h(end+1-(1:nPvb)) = phiB;
    E_h_prev = Q * x_h_prev;
    E_h = Q * x_h;
    rms_e(iter) = sqrt( sum( (E_h - E_h_prev).^2 )) / length(E_h);
    rms_h(iter) = sqrt( sum( (H_h - H_h_prev).^2 )) / length(H_h);
    
    mean_e(iter) = mean(E_h ./ primal.edgeLength);
    mean_h(iter) = mean(H_h ./ dual.edgeLength(1:nPfi));
    
    if mod(iter,10) == 0 || (iter+1) == iterMax
        figure(3)
        plot(timeVec, mean_e, 'DisplayName', '<E>')
        hold on
        plot(timeVec, mean_h, 'DisplayName', '<H>')
        hold off
        ylabel('<E>,<H>')
        legend show
        grid on
        
        figure(1)
        subplot(3,1,1)
        plot(E_h ./ primal.edgeLength,'.')
        grid on
        ylabel('E')
        subplot(3,1,2)
        plot(H_h ./ dual.edgeLength(1:nPfi), '.')
        grid on
        ylabel('H')
        subplot(3,1,3)
        semilogy(timeVec, [rms; rms_e; rms_h])
        ylabel('rms')
        grid on
        
        figure(2)
        clf
        % plot3(primal.vc(:,1), primal.vc(:,2), primal.vc(:,3), '.')
        
        % Plot electric field
        scale = 1;
        pos = primal.edgeCenter;
        vec = E_h ./ primal.edgeLength .* primal.edgeDir;
        idx = 1:nPe;
        quiver3(pos(idx,1), pos(idx,2), pos(idx,3), ...
            vec(idx,1), vec(idx,2), vec(idx,3), scale, 'DisplayName', 'E')
        zlim([-0.1 1.1])
        legend('Location','best')
        %
        % Plot magnetic field
        figure(4)
            pos = primal.faceCenter(1:nPfi,:);
            vec = H_h ./ dual.edgeLength(1:nPfi) .* dual.edgeDir(1:nPfi,:);
%             hold on
            quiver3(pos(:,1), pos(:,2), pos(:,3), ...
                vec(:,1), vec(:,2), vec(:,3), 1, 'DisplayName', 'H')
%             hold off
            grid on
            axis equal
        legend('Location','best')
    end
    
    %     break;
    %     pause
    
    time = time + dt;
    iter = iter + 1;
    
end


%%
return
%% Plot electric field
plot(u, '.')
plot(x_h,'.')
plot(Q*x_h,'.')

%% Reformulated Ampere equation (2nd-order time derivative)
clc
Mnu = spdiags(dual.edgeLength(1:nPf) ./ primal.faceArea, 0, nPf, nPf);
Mnu00 = Mnu(1:nPfi, 1:nPfi);
C00 = primal.f2e(1:nPfi, 1:nPei);
P = [C00 sparse(nPfi,nPvb)];
A = Q' * Meps * Q + P' * Mnu00 * P;
x = zeros(size(A,2),1);
xd = phiT(1);

idxd = false(size(x));
idxd(end+1-(1:nPvt)) = true;
idxf = ~idxd;

Af = A(idxf,idxf);
Ad = A(:,idxd);

rhs = -Ad * xd;
rhsf = rhs(idxf);
xf = Af \ rhsf;

assert(all(size([xf; xd]) == size(x)))
x = [xf; xd];
plot(x,'.')
E_h = Q*x;
vec = E_h ./ primal.edgeLength;
plot(vec, '.')
vec = E_h ./ primal.edgeLength .* primal.edgeDir;
pos = primal.edgeCenter;
quiver3(pos(:,1), pos(:,2), pos(:,3), ...
    vec(:,1), vec(:,2), vec(:,3))


