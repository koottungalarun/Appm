clear
clc

filename = 'Msigma-0.h5';
info = h5info(filename);
% info.Datasets.Name

% Read matrix from file
r = h5read(filename, '/MsigmaInner_rowIdx');
c = h5read(filename, '/MsigmaInner_colIdx');
v = h5read(filename, '/MsigmaInner_values');
Msigma = sparse(r+1,c+1,v);

% Define symmetric matrix, this is numerically more stable
M = 0.5 * (Msigma + Msigma');

issymmetric(Msigma)
issymmetric(M)

% Get difference to original matrix
D = (Msigma - M);
tol = 1e-15 * max(abs(v));
% Show matrix shape
figure(1)
spy(D > tol)

%% Cancel values smaller than a threshold
% tic
[r,c,v] = find(D);
idx = abs(v) < tol;
D(r(idx),c(idx)) = 0;
% D(abs(D) > 0 & abs(D) < tol) = 0;
% toc;

% Get values
[~,~,v] = find(D);
min(abs(v))
figure(2)
histogram(abs(v))

% %% Cut most upper part of matrix
% D = D(1:360,1:360);
% D = triu(D); % get upper triangular part
% spy(D)

% get coefficients
[r,c,v] = find(D);
% plot(unique(sort(r)), '.')

%%
faceIdx = unique(sort(r));
% faceIdx = 1;

primalMesh = readMesh('primal');
f2e = primalMesh.f2e;
e2v = primalMesh.e2v;
f2v = primalMesh.f2v;

vidx = f2v(faceIdx(1),:);

[~,edgeIdx] = find(f2e(faceIdx,:));
% fprintf('edgeIdx: \n')
% disp(edgeIdx)
% for edx = edgeIdx
%     [~,vidx] = find(e2v(edx,:))
% end

vc = primalMesh.vc;

figure(3)
clf
patch('Faces', f2v(faceIdx,:), 'Vertices', vc, 'FaceColor', 'r', 'FaceAlpha', 0.1)
axis equal



