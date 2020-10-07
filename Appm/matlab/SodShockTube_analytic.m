% Analytic solution of Sod shock tube problem
% 
% References:
% 
% Astrophysics; Scaling with speed of sound
% https://www.astro.uu.se/~hoefner/astro/teach/ch10.pdf
% 
% Wikipedia: 
% https://en.wikipedia.org/wiki/Rankine%E2%80%93Hugoniot_conditions
% https://en.wikipedia.org/wiki/Sod_shock_tube
% and
% https://physics.stackexchange.com/questions/423758/how-to-get-exact-solution-to-sod-shock-tube-test
% 
% Dimensional
% https://www3.nd.edu/~gtryggva/CFD-Course2017/Lecture-10-2017.pdf
% 
clear
clc

% Data vectors
tmax = 1;
xR =  2.5;
xL = -xR;
x = linspace(xL, xR, 200); % spatial coordinate
rho = zeros(size(x));      % density
p = zeros(size(x));        % pressure
u = zeros(size(x));        % velocity

% gamma = 7/5;
gamma = 5/3;
Gamma = (gamma-1)/(gamma+1);
beta = (gamma-1)/(2*gamma);

% Left state
% pL = 10/gamma;
% rhoL = 8;
pL = 1;
rhoL = 1;
uL = 0;
cL = sqrt(gamma * pL / rhoL);

% Right state
% pR = 1/gamma;
% rhoR = 1;
pR = 0.1;
rhoR = 0.125;
uR = 0;
cR = sqrt(gamma * pR / rhoR);

x0 = 1/2 * (xL + xR);

pfun1 = @(p) (p - pR) .* sqrt( (1-Gamma) ./ (rhoR * (p + Gamma * pR)) );
pfun2 = @(p) (pL^beta - p.^beta) .* sqrt( ((1 - Gamma^2)*pL.^(1/gamma)) ./ (Gamma^2 * rhoL) );
pfun = @(pp) pfun1(pp) - pfun2(pp);
pInit = 0.5*(pL+pR);
p3 = fzero(pfun,pInit);

rho3 = rhoL * (p3/pL)^(1/gamma);
u3 = uR + (p3 - pR) / sqrt(rhoR/2 * ( (gamma+1) * p3 + (gamma-1)* pR));
u4 = u3;
p4 = p3;

u2   = @(xx,tt) 2/(gamma+1) * (cL + (xx - x0)/tt);
rho2 = @(xx,tt) rhoL * (1 - (gamma-1)/2 * u2(xx,tt)/cL).^(2/(gamma-1));
p2   = @(xx,tt)  pL  * (1 - (gamma-1)/2 * u2(xx,tt)/cL).^(2*gamma/(gamma-1));

rho4 = rhoR * (p4 + Gamma*pR) / (pR + Gamma*p4);


c3 = sqrt(gamma * p3 / rho3);
s_contact = u3;
s_shock = uR + cR * sqrt( ((gamma-1) + (gamma+1) * (p4 / pR)) / (2*gamma) );
s_rfl = -cL;       % speed of  left edge of rarefaction wave
s_rfr = (u4 - c3); % speed of right edge of rarefaction wave

assert(all(diff([s_rfl s_rfr s_contact s_shock]) > 0))

% tvec = [0.01 : 0.01 : tmax];
tvec = tmax;
for t = tvec
% Position of rarefaction wave and discontinuities
x12 = -cL * t;
x23 = s_rfr * t;
x34 = s_contact * t;
x45 = s_shock * t;
xList = [xL x12 x23 x34 x45 xR];
xx = unique(sort([x xList]));
assert(all( diff(xList) > 0 ))


% Left state
x = [];
rho = [];
p = [];
u = [];

idx = xx <= x12;
x = [x xx(idx)];
rho = [rho rhoL*ones(size(xx(idx)))];
p = [p pL*ones(size(xx(idx)))];
u = [u uL*ones(size(xx(idx)))];


% From Left state to rarefaction wave
idx = xx >= x12 & xx <= x23;
x = [x xx(idx)];
rho = [rho rho2(xx(idx),t)];
p = [p p2(xx(idx), t)];
u = [u u2(xx(idx), t)];

% From rarefaction wave to contact discontinuity
idx = xx >= x23 & xx <= x34;
x = [x xx(idx)];
rho = [rho rho3*ones(size(xx(idx)))];
p = [p p3*ones(size(xx(idx)))];
u = [u u3*ones(size(xx(idx)))];

% From contact discontinuity to shock wave
idx = xx >= x34 & xx <= x45;
x = [x xx(idx)];
rho = [rho rho4*ones(size(xx(idx)))];
p = [p p4*ones(size(xx(idx)))];
u = [u u4*ones(size(xx(idx)))];

% % Right state
idx = xx >= x45;
x = [x xx(idx)];
rho = [rho rhoR*ones(size(xx(idx)))];
p = [p pR*ones(size(xx(idx)))];
u = [u uR*ones(size(xx(idx)))];

assert(length(x) == length(rho))
assert(length(x) == length(p))
assert(length(x) == length(p))

figure(1)
subplot(3,1,1)
% density
plot(x,rho)
title('density')

subplot(3,1,2)
% pressure
plot(x,p)
title('pressure')

subplot(3,1,3)
% velocity
plot(x, u)
title('velocity')

end

filename = 'shockTubeAnalytic.mat';
fprintf('Save data to file: %s\n', filename)
save(filename, 'x', 'p', 'rho', 'u', 't', 'gamma')
