function mesh = readMesh(prefix)

% Check if filename exists
filename = strcat(prefix, '-mesh.h5');
assert(exist(filename, 'file') > 0);

% Get list of datasets
info = h5info(filename);
nDatasets = length(info.Datasets);
datasetNames = {};
for i = 1 : nDatasets
    datasetNames{i} = info.Datasets(i).Name;
end

mesh.vc = h5read(filename, '/vertexPos')';
mesh.isBoundaryVertex = h5read(filename, '/isBoundaryVertex');
mesh.isBoundaryVertex = logical(mesh.isBoundaryVertex);
mesh.vertexType = h5read(filename, '/vertexType');

mesh.edgeType = h5read(filename, '/edgeType');
mesh.edgeLength = h5read(filename, '/edgeLength');

mesh.faceArea = h5read(filename, '/faceArea');
mesh.isFaceBoundary = h5read(filename, '/isFaceBoundary')';
mesh.isFaceBoundary = logical(mesh.isFaceBoundary);
mesh.faceCenter = h5read(filename, '/faceCenter')';
mesh.faceNormal = h5read(filename, '/faceNormal')';


mesh.cellCenter = readDataset(filename, 'cellCenter')';
mesh.cellVolume = readDataset(filename, 'cellVolume');
mesh.cellType = readDataset(filename, 'cellFluidType');
% mesh.cellCenter = h5read(filename, '/cellCenter')';
% mesh.cellVolume = h5read(filename, '/cellVolume');
% mesh.cellType = h5read(filename, '/cellFluidType');

% filename = strcat(prefix, '-f2e.dat');
% if exist(filename, 'file')
isMatlabOlderThanR2020a = string(version('-Release')) < "2020a";
    dataname = '/f2e';
    r = h5read(filename, strcat(dataname, '_rowIdx'));
    c = h5read(filename, strcat(dataname, '_colIdx'));
    v = h5read(filename, strcat(dataname, '_values'));
    if isMatlabOlderThanR2020a
        mesh.f2e = sparse(double(r+1),double(c+1),v);
    else
        mesh.f2e = sparse(r+1,c+1,v);
    end
% end

% filename = strcat(prefix, '-e2v.dat');
% if exist(filename, 'file')
    dataname = '/e2v';
    r = h5read(filename, strcat(dataname, '_rowIdx'));
    c = h5read(filename, strcat(dataname, '_colIdx'));
    v = h5read(filename, strcat(dataname, '_values'));
    if isMatlabOlderThanR2020a
        mesh.e2v = sparse(double(r+1),double(c+1),v);
    else
        mesh.e2v = sparse(r+1,c+1,v);
    end
% end

% filename = strcat(prefix, '-c2f.dat');
% if exist(filename, 'file')
    dataname = '/c2f';
    r = h5read(filename, strcat(dataname, '_rowIdx'));
    c = h5read(filename, strcat(dataname, '_colIdx'));
    v = h5read(filename, strcat(dataname, '_values'));
    if isMatlabOlderThanR2020a
        mesh.c2f = sparse(double(r+1),double(c+1),v);
    else
        mesh.c2f = sparse(r+1,c+1,v);
    end
% end

filename = strcat(prefix, '-f2v.dat');
mesh.f2v = readIfFileExists(filename);
if ~isempty(mesh.f2v)
    mesh.f2v(mesh.f2v(:) < 0) = nan;
    mesh.f2v = mesh.f2v + 1;
end



% edge direction
if isfield(mesh, 'e2v') && isfield(mesh, 'vc')
    [r,c] = find(mesh.e2v == -1);
    idxA(r) = c;
    [r,c] = find(mesh.e2v == 1);
    idxB(r) = c;
    mesh.edgeDir = mesh.vc(idxB,:) - mesh.vc(idxA,:);
    mesh.edgePos = mesh.vc(idxA,:);
    mesh.edgeCenter = 0.5 * (mesh.vc(idxA, :) + mesh.vc(idxB,:));
end

%% Nested functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    function data = readDataset(filename, datasetname)
        % Check if dataset exists
        data = [];
        if any(strcmp(datasetNames, datasetname))
            % Read data
            data = h5read(filename, strcat('/', datasetname));
        else
%             error('Dataset name %s does not exist in file %s.', datasetname, filename)
        end
    end



    function data = readIfFileExists(filename)
        
        if exist(filename, 'file')
            data = importdata(filename);
        else
            data = [];
        end
        
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
