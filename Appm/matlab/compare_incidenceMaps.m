clear
primal = readMesh('primal');
dual = readMesh('dual');

%% Compare face-to-edge maps

figure(1)
subplot(1,2,1)
imagesc(primal.f2e), colorbar

nPe = size(primal.f2e, 2);
nPf = size(primal.f2e, 1);

subplot(1,2,2)
imagesc(dual.f2e(1:nPe,1:nPf)'), colorbar

shg

T = primal.f2e - dual.f2e(1:nPe,1:nPf)';
nnzT = nnz(T);
assert(nnz(T) == 0)

%% Compare primal edge-to-vertexmap with dual cell-to-face map
clf
subplot(1,2,1)
imagesc(primal.e2v), colorbar

nPv = size(primal.e2v,2);
nPc = size(primal.c2f,1);

subplot(1,2,2)
imagesc(-dual.c2f(1:nPv, 1:nPe)'), colorbar
shg
