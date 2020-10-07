% Check if the matrix Msigma is symmetric positive definite


%% Matrix with form factors of unit normals
M_rjni = readSparseMatrix('M_rjni.dat');
[~,~,v] = find(M_rjni);
fprintf('min,max(M_rjni): %e, %e \n', min(abs(v)), max(abs(v))); % analyse size of values in T
vsorted = sort(abs(v));
histogram(log10(vsorted))

assert(min(abs(v)) > 1e-2, 'Form factors rj.ni is bad')


%% Load data

Q = readSparseMatrix('Q.dat');
spy(Q)

Msigma = readSparseMatrix('Msigma.dat');
spy(Msigma)

[r,c,v] = find(Msigma);
vsorted = sort(abs(v));
histogram(log10(vsorted));

% Is Msigma symmetric, approximately?
fprintf('Msigma symmetric: %d \n', issymmetric(Msigma));

T = Msigma - Msigma';
fprintf('nnz(T): %d \n', nnz(T))
[r,c,v] = find(T); % get non-zero entries in T
fprintf('min(T),max(T): %e, %e \n', min(abs(v)), max(abs(v))); % analyse size of values in T
assert(max(abs(v)) < 1e-14);

fprintf('Msigma is symmetric up to error = %e \n', max(abs(v)));

%% Make Msigma symmetric

Z = 0.5 * (Msigma + Msigma'); % make symmetric
result(1) = issymmetric(Z);

M = Q' * (0.5 * (Msigma + Msigma')) * Q; % not symmetric
result(2) = issymmetric(M);

M = 0.5 * Q' * (Msigma * Q) + 0.5 * (Q' * Msigma') * Q; % is symmetric!
result(3) = issymmetric(M);

disp(result)

primalMesh = readMesh('primal');
nEdges = size(Q,1);
nDof = size(Q,2);
nTerminalVertices = 218;
nF = nDof - nTerminalVertices;

idxF = [
    true(nF,1);
    false(nTerminalVertices,1)
    ];
idxD = ~idxF;

%     Af = M(idxF,idxF);

%%
for n = 10:55
    disp(n)
    A = M(1:n,1:n);
    lambda = eig(full(A),eye(size(A)),'qz');
%     lambda(1:10)'
    min(lambda)
    idx = lambda < 0;
    v = 1 : length(lambda);
    semilogy(v, abs(lambda), v(idx), abs(lambda(idx)), '*')
    grid on
%     chol(A);
    pause
end



return
T = Z - Msigma; % difference to actual matrix (check for numerical inaccuracies)
[r,c,v] = find(T); % get values
nnz(T) / nnz(Msigma)
vsorted = sort(abs(v));
histogram(log10(vsorted))
%%

M = Q' * Z * Q;
issymmetric(M)


%% Check if M is approximately symmetric

T = M - M'; % T is zero if M is symmetric
fprintf('nnz(T): %d \n', nnz(T))
[r,c,v] = find(T); % get non-zero entries in T
fprintf('min(T),max(T): %e, %e \n', min(abs(v)), max(abs(v))); % analyse size of values in T

tol = 5e-15;
spy(abs(T) > 4e-15)


% Edit data in matrix M
figure(1)
[r,c,v] = find(M);
vsorted = sort(abs(v));
histogram(log10(vsorted))

figure(3)
spy(abs(M) < 1e-10 & M ~= 0)

figure(2)
M(abs(M) < 1e-10) = 0;
[r,c,v] = find(M);
vsorted = sort(abs(v));
histogram(log10(vsorted))

