% Plot mesh

coords = importdata('primal-vertices.dat');
e2v = readSparseMatrix('primal-e2v.dat');

[row,col] = find(e2v == -1);
[~,idx] = sort(row);
idxA = col(idx);
[row,col] = find(e2v == +1);
[~,idx] = sort(row);
idxB = col(idx);

A = coords(:,idxA);
B = coords(:,idxB);
edgeDir = B - A;

quiver3(A(1,:), A(2,:), A(3,:), ...
    edgeDir(1,:), edgeDir(2,:), edgeDir(3,:), 0)


