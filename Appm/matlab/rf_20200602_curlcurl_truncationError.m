clear
clc

primal = readMesh('primal');
M = readSparseMatrix('M.dat');
Q = readSparseMatrix('Q.dat');


nDof = size(M,2);
nD = sum(primal.vertexType == 2);
nFree = nDof - nD;

Mf = M(1:nFree, 1:nFree);
Md = M(:, (end+1 - nD) : end);
assert(size(Md,2) == nD);

phiA = 0 * ones(nD/2,1);
phiB = 50 * ones(nD/2,1);
xD = [phiA; 
    phiB];

rhs = zeros(nDof,1);
rhs = rhs - Md * xD;

rhsFree = rhs(1:nFree);
xf = Mf \ rhsFree; % mldivide allows for sparse matrix
% xf = linsolve(full(Mf),rhsFree); % linsolve requires a full matrix
x = [xf; xD];

E_h = Q * x;

Emax = max(abs(E_h));
%% cleaning of low-order bits
scale = 1e7;
temp = Emax * scale;
EE = E_h + temp;
EE = EE - temp;

temp = [E_h EE];

edgeDir = primal.edgeDir;
isEdgeZparallel = false(size(E_h));
nPe = size(primal.e2v,1);
isEdgeZparallel = edgeDir(:,3) ~= 0;

figure(1)
plot(temp(isEdgeZparallel,:))

figure(2)
plot(temp(~isEdgeZparallel,:))
if all(EE(~isEdgeZparallel) == 0)
    disp('zparallels are zero')
else
    disp('zparallels are NOT zero')
end


