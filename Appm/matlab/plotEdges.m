function plotEdges(mesh, showEdgeLabels, color)

gcf;

% plot edges
[r,c] = find(mesh.e2v == -1);
idxA(r) = c;
[r,c] = find(mesh.e2v == 1);
idxB(r) = c;
edgeDir = mesh.vc(idxB,:) - mesh.vc(idxA,:);
edgeCenter = 1/2 * (mesh.vc(idxA,:) + mesh.vc(idxB,:));
hold on
quiver3(mesh.vc(idxA,1), mesh.vc(idxA,2), mesh.vc(idxA,3), ...
    edgeDir(:,1), edgeDir(:,2), edgeDir(:,3), 0)
if showEdgeLabels
    text(edgeCenter(:,1), edgeCenter(:,2), edgeCenter(:,3), num2cell( 1:size(edgeCenter,1) ), 'Color', color)
end
hold off



end